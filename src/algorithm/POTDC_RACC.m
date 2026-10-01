function result = POTDC_RACC(HB_ctrl, HD_ctrl, para)
%POTDC_RACC POTDC solution of the ATF/correlation hybrid robust ACC model.
%
% The implemented source-paper problem is
%
%   minimize_w   w^H (H_D^H H_D + gamma_D I) w
%   subject to   ||H_B w||_2 - eta_B ||w||_2 >= 1,                 (1)
%
% where eta_B is an ATF-level bright-zone radius and gamma_D is a
% correlation-level dark-zone radius. In the controlled comparison, the
% common ATF radii are mapped as
%
%   eta_B(f)   = nu ||H_B(f)||_F,
%   gamma_D(f) = eta_D(f)^2,  eta_D(f) = nu ||H_D(f)||_F.          (2)
%
% Thus POTDC uses the same normalized uncertainty level as the other robust
% algorithms while retaining the hybrid uncertainty model of the original
% method. Historical callers without explicit common radii keep the former
% gamma_D=10 and eta_B=1e-3*sqrt(trace(H_B^H H_B)) settings.
%
% Optional para.potdc fields:
%   iterations       maximum POTDC iterations (default 20)
%   tolerance        relative objective tolerance (default 1e-8)
%   solver           CVX solver; empty keeps the active CVX solver
%   mosek_threads    MOSEK threads per worker (default 2)
%   parallel_workers frequency-level MATLAB workers (default 1)

    if ndims(HB_ctrl) == 2
        HB_ctrl = reshape(HB_ctrl, size(HB_ctrl, 1), size(HB_ctrl, 2), 1);
    end
    if ndims(HD_ctrl) == 2
        HD_ctrl = reshape(HD_ctrl, size(HD_ctrl, 1), size(HD_ctrl, 2), 1);
    end
    [~, L, num_frequencies] = size(HB_ctrl);
    if size(HD_ctrl, 2) ~= L || size(HD_ctrl, 3) ~= num_frequencies
        error('POTDC:DimensionMismatch', ...
            'HB_ctrl and HD_ctrl must have compatible dimensions.');
    end

    opts = local_options(para);
    use_common_radii = local_has_common_radii(para);
    frequency_results = cell(1, num_frequencies);

    if opts.parallel_workers > 1 && num_frequencies > 1
        local_ensure_parallel_pool(opts.parallel_workers);
        parfor (i = 1:num_frequencies, opts.parallel_workers)
            frequency_results{i} = local_solve_frequency( ...
                HB_ctrl(:, :, i), HD_ctrl(:, :, i), para, opts, ...
                use_common_radii, i);
        end
    else
        for i = 1:num_frequencies
            frequency_results{i} = local_solve_frequency( ...
                HB_ctrl(:, :, i), HD_ctrl(:, :, i), para, opts, ...
                use_common_radii, i);
        end
    end

    w = zeros(L, num_frequencies);
    status = repmat({''}, 1, num_frequencies);
    etaB = nan(1, num_frequencies);
    etaD = nan(1, num_frequencies);
    gammaD = nan(1, num_frequencies);
    rank_fraction = nan(1, num_frequencies);
    rank_ratio = nan(1, num_frequencies);
    iteration_count = zeros(1, num_frequencies);
    objective = nan(1, num_frequencies);
    solve_time_seconds = nan(1, num_frequencies);
    for i = 1:num_frequencies
        item = frequency_results{i};
        w(:, i) = item.w;
        status{i} = item.status;
        etaB(i) = item.etaB;
        etaD(i) = item.etaD;
        gammaD(i) = item.gammaD;
        rank_fraction(i) = item.rank_fraction;
        rank_ratio(i) = item.rank_ratio;
        iteration_count(i) = item.iterations;
        objective(i) = item.objective;
        solve_time_seconds(i) = item.solve_time_seconds;
    end

    result = struct();
    result.w = w;
    result.scale = 1;
    result.status = status;
    result.eta = struct('B', etaB, 'D', etaD);
    result.gamma = struct('D', gammaD);
    result.uses_common_atf_radii = use_common_radii;
    result.uncertainty_mapping = ...
        'POTDC hybrid model: bright eta_B; dark gamma_D=eta_D^2';
    result.rank_fraction = rank_fraction;
    result.rank_ratio = rank_ratio;
    result.iterations = iteration_count;
    result.objective = objective;
    result.solve_time_seconds = solve_time_seconds;
    result.options = opts;
end

function item = local_solve_frequency(HB, HD, para, opts, use_common, i)
    RB = (HB' * HB + (HB' * HB)') / 2;
    if use_common
        etaB = local_resolve_eta(para, 'B', i, HB);
        etaD = local_resolve_eta(para, 'D', i, HD);
        gammaD = etaD^2;
    else
        etaB = 1e-3 * sqrt(max(real(trace(RB)), 0));
        etaD = NaN;
        gammaD = 10;
    end

    % Scale H and both uncertainty quantities together. Solving (1) with
    % H/data_scale changes only the magnitude of w; dividing the returned
    % vector by data_scale restores the original constraint normalization.
    data_scale = max(norm([HB; HD], 'fro'), realmin);
    HBs = HB / data_scale;
    HDs = HD / data_scale;
    etaBs = etaB / data_scale;
    gammaDs = gammaD / data_scale^2;

    timer = tic;
    [w_scaled, diagnostics] = local_potdc_solver( ...
        HBs, HDs, gammaDs, etaBs, opts);
    diagnostics.solve_time_seconds = toc(timer);
    item = diagnostics;
    item.w = w_scaled / data_scale;
    item.etaB = etaB;
    item.etaD = etaD;
    item.gammaD = gammaD;

    fprintf(['POTDC frequency %d: %s, iterations %d, rank fraction ', ...
        '%.4f, %.1f s.\n'], i, item.status, item.iterations, ...
        item.rank_fraction, item.solve_time_seconds);
end

function [w_opt, diagnostics] = local_potdc_solver( ...
        HB, HD, gammaD, etaB, opts)
    L = size(HB, 2);
    RB = (HB' * HB + (HB' * HB)') / 2;
    RD = (HD' * HD + (HD' * HD)') / 2;
    [V, D] = eig(RB);
    [max_eigenvalue_RB, idx] = max(real(diag(D)));
    bright_margin = sqrt(max(max_eigenvalue_RB, 0)) - etaB;
    if bright_margin <= 0
        error('POTDC:EmptyBrightRobustSet', ...
            ['eta_B=%.4g is not smaller than sqrt(lambda_max(R_B))=', ...
             '%.4g; the POTDC initialization is infeasible.'], ...
            etaB, sqrt(max(max_eigenvalue_RB, 0)));
    end

    dark_matrix = RD + gammaD * eye(L);
    numerical_floor = max(1e-12 * max(real(trace(dark_matrix)) / L, 1), eps);
    dark_matrix = dark_matrix + numerical_floor * eye(L);
    generalized_values = eig(RB, dark_matrix);
    max_eigenvalue_RBD = max(real(generalized_values));

    alpha_lower = 1 / (1 - etaB / sqrt(max_eigenvalue_RB))^2;
    w0 = (1.001 / bright_margin) * V(:, idx);
    alpha_upper = max_eigenvalue_RBD * real(w0' * dark_matrix * w0);
    alpha_upper = max(alpha_upper, alpha_lower * (1 + 1e-8));
    alpha_center = (alpha_upper + alpha_lower) / 2;

    W_current = [];
    objective_previous = Inf;
    objective_current = NaN;
    final_status = 'Not run';
    completed_iterations = 0;
    for iteration = 1:opts.iterations
        cvx_begin sdp quiet
            if ~isempty(opts.solver)
                cvx_solver(opts.solver)
                if startsWith(lower(opts.solver), 'mosek') && ...
                        ~isempty(opts.mosek_threads)
                    cvx_solver_settings( ...
                        'MSK_IPAR_NUM_THREADS', opts.mosek_threads)
                end
            end
            variable W(L, L) hermitian semidefinite
            variable alpha_opt nonnegative
            minimize(real(trace(dark_matrix * W)))
            subject to
                real(trace(RB * W)) == alpha_opt;
                alpha_lower <= alpha_opt;
                alpha_opt <= alpha_upper;
                etaB^2 * real(trace(W)) + (sqrt(alpha_center) - 1) + ...
                    alpha_opt * (1 / sqrt(alpha_center) - 1) <= 0;
        cvx_end

        final_status = cvx_status;
        completed_iterations = iteration;
        if ~local_status_is_solved(final_status)
            break;
        end
        W_current = full((W + W') / 2);
        alpha_center = full(alpha_opt);
        objective_current = real(trace(dark_matrix * W_current));
        relative_change = abs(objective_previous - objective_current) / ...
            max(abs(objective_previous), 1);
        if iteration >= 2 && relative_change <= opts.tolerance
            break;
        end
        objective_previous = objective_current;
    end

    if isempty(W_current)
        w_opt = zeros(L, 1);
        rank_fraction = NaN;
        rank_ratio = NaN;
    else
        [eigenvectors, eigenvalue_matrix] = eig(W_current);
        eigenvalues = max(real(diag(eigenvalue_matrix)), 0);
        [lambda_maximum, maximum_index] = max(eigenvalues);
        sorted_values = sort(eigenvalues, 'descend');
        rank_fraction = lambda_maximum / max(sum(eigenvalues), realmin);
        if numel(sorted_values) > 1
            rank_ratio = sorted_values(1) / max(sorted_values(2), eps);
        else
            rank_ratio = Inf;
        end
        w_opt = sqrt(lambda_maximum) * eigenvectors(:, maximum_index);
    end

    diagnostics = struct('status', final_status, ...
        'iterations', completed_iterations, 'objective', objective_current, ...
        'rank_fraction', rank_fraction, 'rank_ratio', rank_ratio, ...
        'solve_time_seconds', NaN);
end

function tf = local_has_common_radii(para)
    tf = (isfield(para, 'epsilon') && ...
        isfield(para.epsilon, 'B') && isfield(para.epsilon, 'D')) || ...
        (isfield(para, 'full_racc') && ...
        (isfield(para.full_racc, 'eta') || ...
         isfield(para.full_racc, 'relative_eta')));
end

function eta = local_resolve_eta(para, zone, frequency_index, H)
    if isfield(para, 'full_racc') && isfield(para.full_racc, 'eta') && ...
            isfield(para.full_racc.eta, zone)
        eta = local_pick_value(para.full_racc.eta.(zone), frequency_index);
    elseif isfield(para, 'epsilon') && isfield(para.epsilon, zone)
        eta = local_pick_value(para.epsilon.(zone), frequency_index);
    elseif isfield(para, 'full_racc') && ...
            isfield(para.full_racc, 'relative_eta')
        relative_source = para.full_racc.relative_eta;
        if isstruct(relative_source)
            relative_source = relative_source.(zone);
        end
        eta = local_pick_value(relative_source, frequency_index) * ...
            norm(H, 'fro');
    else
        error('POTDC:MissingUncertaintyRadius', ...
            'The common ATF uncertainty radius is incomplete.');
    end
    validateattributes(eta, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
end

function value = local_pick_value(source, frequency_index)
    if isscalar(source)
        value = source;
    elseif numel(source) >= frequency_index
        value = source(frequency_index);
    else
        error('POTDC:UncertaintyLengthMismatch', ...
            'No uncertainty radius is available for frequency index %d.', ...
            frequency_index);
    end
end

function opts = local_options(para)
    opts = struct('iterations', 20, 'tolerance', 1e-8, 'solver', '', ...
        'mosek_threads', 2, 'parallel_workers', 1);
    if isfield(para, 'potdc')
        names = fieldnames(opts);
        for k = 1:numel(names)
            if isfield(para.potdc, names{k})
                opts.(names{k}) = para.potdc.(names{k});
            end
        end
    end
    opts.solver = char(string(opts.solver));
    validateattributes(opts.iterations, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(opts.tolerance, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(opts.parallel_workers, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
end

function local_ensure_parallel_pool(requested_workers)
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', requested_workers);
    elseif pool.NumWorkers < requested_workers
        warning('POTDC:ExistingPoolSmallerThanRequested', ...
            ['The existing pool has %d workers; requested %d. The existing ', ...
             'pool is retained.'], pool.NumWorkers, requested_workers);
    end
end

function solved = local_status_is_solved(status)
    solved = strcmp(status, 'Solved') || strcmp(status, 'Inaccurate/Solved');
end
