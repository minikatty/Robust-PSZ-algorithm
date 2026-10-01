function result = RACC_PM_GlobalMatched(HB_ctrl, HD_ctrl, H_desired, para, fre_indices)
%RACC_PM_GLOBALMATCHED Raw, un-decomposed counterpart of RACC_PM_Sub.
%   This implementation is intended only for the controlled Schur-
%   decomposition timing/consistency ablation.  It uses the same ATF-level
%   radii, normalization retries, objective, power bound, solver, and rank
%   recovery rule as the maintained no-suffix RACC_PM_Sub implementation.
%
%   The raw robust LMIs have orders MB*L+1 and MD*L+1.  At the full
%   96-point, 48-loudspeaker configuration, each order is 4609 and CVX can
%   require prohibitive memory.  A configurable guard prevents accidental
%   allocation of that problem unless explicitly overridden.

    [MB, L, num_frequency_points] = size(HB_ctrl);
    MD = size(HD_ctrl, 1);
    if nargin < 5 || isempty(fre_indices)
        loop_target = 1:num_frequency_points;
    else
        loop_target = fre_indices;
    end

    rho = para.rho;
    mu = para.mu;
    gamma = mu + rho .* para.alpha;
    epsilon_B = para.scale .* para.epsilon.B;
    epsilon_D = para.scale .* para.epsilon.D;
    power_bound = 1e4;
    normalization_multipliers = [0.1, 1, 5, 10, 0.01];
    mosek_threads = [];
    rank_fraction_threshold = 0.98;
    max_raw_lmi_order = 2000;
    allow_large_raw_lmi = false;
    if isfield(para, 'racc_pm')
        options = para.racc_pm;
        if isfield(options, 'mosek_threads')
            mosek_threads = options.mosek_threads;
        end
        if isfield(options, 'rank_fraction_threshold')
            rank_fraction_threshold = options.rank_fraction_threshold;
        end
    end
    if isfield(para, 'racc_pm_global')
        options = para.racc_pm_global;
        if isfield(options, 'max_raw_lmi_order')
            max_raw_lmi_order = options.max_raw_lmi_order;
        end
        if isfield(options, 'allow_large_raw_lmi')
            allow_large_raw_lmi = logical(options.allow_large_raw_lmi);
        end
        if isfield(options, 'normalization_multipliers')
            normalization_multipliers = options.normalization_multipliers;
        end
    end

    largest_raw_lmi_order = max(MB, MD) * L + 1;
    if largest_raw_lmi_order > max_raw_lmi_order && ~allow_large_raw_lmi
        error('RACC_PM_GlobalMatched:RawProblemTooLarge', ...
            ['Raw LMI order %d exceeds the safety limit %d. The full ', ...
             'problem is intentionally not allocated because it can exhaust ', ...
             'memory. Set para.racc_pm_global.allow_large_raw_lmi=true ', ...
             'only for an explicitly monitored run.'], ...
            largest_raw_lmi_order, max_raw_lmi_order);
    end

    weights = zeros(L, num_frequency_points);
    status = strings(1, num_frequency_points);
    objective_value = nan(1, num_frequency_points);
    rank_fraction = nan(1, num_frequency_points);
    selected_multiplier = nan(1, num_frequency_points);
    attempt_status = cell(1, num_frequency_points);
    cvx_clear;

    for frequency_index = loop_target(:).'
        HB_raw = HB_ctrl(:, :, frequency_index);
        HD_raw = HD_ctrl(:, :, frequency_index);
        desired_raw = H_desired(:, frequency_index);
        eta_B = local_pick_value(epsilon_B, frequency_index);
        eta_D = local_pick_value(epsilon_D, frequency_index);
        gamma_i = local_pick_value(gamma, frequency_index);

        base_scale = norm([HB_raw; HD_raw], 'fro') * sqrt(rho);
        if base_scale < 1e-10
            base_scale = 1;
        end
        statuses = strings(1, numel(normalization_multipliers));
        for retry_index = 1:numel(normalization_multipliers)
            multiplier = normalization_multipliers(retry_index);
            current_scale = base_scale * multiplier;
            H_i = struct('B', HB_raw / current_scale, ...
                'D', HD_raw / current_scale);
            desired_i = desired_raw / current_scale;
            eta_i = struct('B', eta_B / current_scale, ...
                'D', eta_D / current_scale);
            solver_options = struct('rho', rho, 'gamma', gamma_i, ...
                'power_bound', power_bound, 'L', L, 'MB', MB, 'MD', MD, ...
                'mosek_threads', mosek_threads, ...
                'rank_fraction_threshold', rank_fraction_threshold);
            [candidate, candidate_status, candidate_objective, ...
                candidate_rank_fraction] = local_solve_raw_lmi( ...
                H_i, desired_i, eta_i, solver_options);
            statuses(retry_index) = string(candidate_status);
            if local_is_solved(candidate_status)
                weights(:, frequency_index) = candidate;
                status(frequency_index) = string(candidate_status);
                objective_value(frequency_index) = candidate_objective;
                rank_fraction(frequency_index) = candidate_rank_fraction;
                selected_multiplier(frequency_index) = multiplier;
                break;
            end
        end
        attempt_status{frequency_index} = statuses;
        if strlength(status(frequency_index)) == 0
            status(frequency_index) = statuses(end);
        end
    end

    result = struct('w', weights, 'status', status, ...
        'objective_value', objective_value, ...
        'rank_fraction', rank_fraction, ...
        'normalization_multiplier', selected_multiplier, ...
        'normalization_attempt_status', {attempt_status}, ...
        'raw_lmi_order_B', MB * L + 1, ...
        'raw_lmi_order_D', MD * L + 1, ...
        'formulation', 'global-raw');
end

function value = local_pick_value(source, frequency_index)
    if isscalar(source)
        value = source;
    else
        value = source(frequency_index);
    end
end

function solved = local_is_solved(status)
    solved = strcmp(status, 'Solved') || ...
        strcmp(status, 'Inaccurate/Solved');
end

function [recovered_w, status, objective_value, rank_fraction] = ...
        local_solve_raw_lmi(H_i, desired, eta, options)
    rho = options.rho;
    gamma = options.gamma;
    L = options.L;
    MB = options.MB;
    MD = options.MD;
    power_bound = options.power_bound;
    noise_floor = 1e-8;

    cvx_begin sdp quiet
        cvx_solver mosek
        if ~isempty(options.mosek_threads)
            cvx_solver_settings('MSK_IPAR_NUM_THREADS', ...
                options.mosek_threads);
        end
        variable W_tilde(L + 1, L + 1) hermitian
        variable tB
        variable tD
        variable tau_B nonnegative
        variable tau_D nonnegative
        expression W(L, L)
        expression w(L, 1)
        expression Q_B(MB * L, MB * L)
        expression Q_D(MD * L, MD * L)
        expression u_B(MB * L, 1)
        expression u_D(MD * L, 1)

        W = W_tilde(1:L, 1:L);
        w = W_tilde(1:L, L + 1);
        Q_B = (1 - rho) * kron(conj(W), eye(MB));
        Q_D = gamma * kron(conj(W), eye(MD));
        u_B = reshape((1 - rho) * H_i.B * W - desired * w', ...
            MB * L, 1);
        u_D = reshape(gamma * H_i.D * W, MD * L, 1);
        cB = (1 - rho) * real(trace(H_i.B * W * H_i.B')) ...
            - 2 * real(desired' * H_i.B * w) ...
            + desired' * desired + noise_floor - tB;
        cD = gamma * real(trace(H_i.D * W * H_i.D')) - tD;

        minimize(tB + tD)
        subject to
            W_tilde >= 0;
            W_tilde(L + 1, L + 1) == 1;
            real(trace(W)) <= power_bound;
            [tau_B * eye(MB * L) - Q_B, -u_B; ...
             -u_B', -cB - tau_B * eta.B^2] >= 0;
            [tau_D * eye(MD * L) - Q_D, -u_D; ...
             -u_D', -cD - tau_D * eta.D^2] >= 0;
    cvx_end

    status = cvx_status;
    objective_value = cvx_optval;
    rank_fraction = NaN;
    recovered_w = zeros(L, 1);
    if ~local_is_solved(status)
        return;
    end

    W_tilde = (W_tilde + W_tilde') / 2;
    [eigenvectors, eigenvalues] = eig(full(W_tilde), 'vector');
    eigenvalues = max(real(eigenvalues), 0);
    [eigenvalues, order] = sort(eigenvalues, 'descend');
    eigenvectors = eigenvectors(:, order);
    rank_fraction = eigenvalues(1) / max(sum(eigenvalues), realmin);
    lifted_mean = W_tilde(1:L, L + 1);
    if abs(eigenvectors(L + 1, 1)) > 1e-10
        principal_w = eigenvectors(1:L, 1) / eigenvectors(L + 1, 1);
    else
        principal_w = lifted_mean;
    end
    if norm(principal_w)^2 > power_bound
        principal_w = principal_w * sqrt(power_bound) / norm(principal_w);
    end
    if rank_fraction >= options.rank_fraction_threshold
        recovered_w = principal_w;
        return;
    end

    covariance = W_tilde(1:L, 1:L) - lifted_mean * lifted_mean';
    covariance = (covariance + covariance') / 2;
    [basis, covariance_eigenvalues] = eig(full(covariance), 'vector');
    covariance_root = basis * diag(sqrt(max(real(covariance_eigenvalues), 0)));
    recovered_w = principal_w;
    best_objective = local_nominal_objective( ...
        recovered_w, H_i, desired, rho, gamma);
    for randomization_index = 1:1000
        random_vector = (randn(L, 1) + 1i * randn(L, 1)) / sqrt(2);
        candidate = lifted_mean + covariance_root * random_vector;
        if norm(candidate)^2 > power_bound
            candidate = candidate * sqrt(power_bound) / norm(candidate);
        end
        candidate_objective = local_nominal_objective( ...
            candidate, H_i, desired, rho, gamma);
        if candidate_objective < best_objective
            recovered_w = candidate;
            best_objective = candidate_objective;
        end
    end
end

function value = local_nominal_objective(w, H_i, desired, rho, gamma)
    bright_pressure = H_i.B * w;
    bright_error = desired - bright_pressure;
    dark_pressure = H_i.D * w;
    value = real(bright_error' * bright_error ...
        - rho * (bright_pressure' * bright_pressure) ...
        + gamma * (dark_pressure' * dark_pressure));
end
