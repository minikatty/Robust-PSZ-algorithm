function res = FullCrossTermRACC(HB_ctrl, HD_ctrl, para, fre_indices)
%FULLCROSSTERMRACC Robust ACC retaining all ATF-induced cross terms.
%
% For W=ww^H, E_Z(W,DeltaH_Z) denotes the loaded zone energy
%
%   E_Z = trace((H_Z+DeltaH_Z)*W*(H_Z+DeltaH_Z)^H)
%         + delta_Z*trace(W),   ||DeltaH_Z||_F <= eta_Z.
%
% The rank-constrained robust ACC problem is
%
%   maximize_{W,s_B,s_D}  (M_D/M_B) s_B/s_D
%   subject to            min_DeltaH_B E_B >= s_B,
%                         max_DeltaH_D E_D <= s_D,
%                         W >= 0, rank(W)=1.
%
% Dropping rank(W)=1 gives the SDP relaxation. With
%
%   Q_Z(W)=W^T kron I_M,  u_Z(W)=vec(H_Z W),
%
% the complex S-procedure produces the exact raw robust LMIs. The
% decomposed form used in the experiment is their Schur-equivalent row-block
% representation: for each bright row m and dark row k,
%
%   [ W^T+lambda_B I,  (h_B,m W)^T ] >= 0,
%   [       *,               d_B,m ]
%   sum_m d_B,m <= c_B-s_B-lambda_B eta_B^2,
%
%   [ lambda_D I-W^T, -(h_D,k W)^T ] >= 0,
%   [          *,             d_D,k ]
%   sum_k d_D,k <= s_D-c_D-lambda_D eta_D^2,
%
% where c_Z=trace(H_Z W H_Z^H)+delta_Z trace(W). MATLAB W.' is the
% ordinary transpose required by W^T; W' is the conjugate transpose.
% Two equivalent solution methods are available:
%
%   'normalized' (default): fix s_D to a positive data-scaled reference by
%                           homogeneity and maximize s_B; this needs only
%                           one SDP per frequency.
%   'bisection': retain trace(W)=1 and bisect over a fixed alpha; this is
%                retained as an independent implementation check.
%
% Two equivalent LMI formulations are also available:
%
%   'decomposed' (default): 2*M LMIs of size L+1; intended for experiments.
%   'raw': one (M*L+1)-LMI per zone; intended only for small validation.
%
% Required uncertainty input (first available form is used):
%   para.full_racc.eta.B/D       absolute Frobenius radii;
%   para.epsilon.B/D             compatibility alias;
%   para.full_racc.relative_eta  scalar or struct, eta_Z=r_Z*||H_Z||_F.
%
% Optional para.full_racc fields:
%   solution_method          'normalized' (default) or 'bisection'
%   formulation             'decomposed' (default) or 'raw'
%   solver                  CVX solver name (decomposed form requires MOSEK)
%   mosek_threads           MOSEK threads used by each MATLAB worker (default 4)
%   bisection_tolerance_db  default 0.10 dB
%   max_bisection_iterations default 14
%   alpha_floor_db          default -80 dB
%   alpha_upper_db          explicit initial upper bound
%   upper_margin_db         default 10 dB over regularized nominal ACC
%   maximum_upper_db        default 160 dB
%   randomization_trials    default 300
%   random_seed             default 20260820
%   recovery_method         'paper' (default) or 'enhanced'
%   refinement_iterations   default 300
%   refinement_seeds        default 8
%   parallel_workers        frequency-level MATLAB workers (default 1)
%   checkpoint_dir          optional per-frequency checkpoint directory
%   resume_from_checkpoints reuse matching checkpoints (default true)
%   normalization_sD        positive fixed s_D; default uses the nominal
%                           mean dark-zone eigenvalue after data scaling
%   normalization_multipliers numerical ATF scaling trials; multiplying H,
%                           eta, and loading consistently leaves the robust
%                           problem unchanged (default 1)
%   retry_unsolved_checkpoints retry matching checkpoints whose status is
%                           not solved (default true)
%   fallback_to_bisection   after all normalization trials fail, use the
%                           trace-normalized fixed-alpha method (default false)
%   diagonal_loading_ratio common loading ratio tau; the regularized zone
%                           energy adds tau*lambda_max(R_Z)*||w||_2^2
%                           (default 0 for backward compatibility)
%   raw_max_dimension       default 800

    if ndims(HB_ctrl) == 2
        HB_ctrl = reshape(HB_ctrl, size(HB_ctrl, 1), size(HB_ctrl, 2), 1);
    end
    if ndims(HD_ctrl) == 2
        HD_ctrl = reshape(HD_ctrl, size(HD_ctrl, 1), size(HD_ctrl, 2), 1);
    end

    [MB, L, num_freqs] = size(HB_ctrl);
    [MD, LD, num_freqs_D] = size(HD_ctrl);
    if LD ~= L || num_freqs_D ~= num_freqs
        error('FullRACC:DimensionMismatch', ...
            'HB_ctrl and HD_ctrl must have compatible loudspeaker/frequency dimensions.');
    end
    if nargin < 4 || isempty(fre_indices)
        fre_indices = 1:num_freqs;
    end
    if any(fre_indices < 1) || any(fre_indices > num_freqs) || ...
            any(fre_indices ~= round(fre_indices))
        error('FullRACC:InvalidFrequencyIndices', ...
            'fre_indices must contain valid integer indices into the third dimension.');
    end

    opts = local_options(para);
    if strcmpi(opts.formulation, 'decomposed') && ...
            ~startsWith(lower(opts.solver), 'mosek')
        error('FullRACC:DecomposedRequiresMosek', ...
            ['The decomposed Full-RACC LMI must be solved with MOSEK in ', ...
             'this project; SDPT3/SeDuMi may exhaust memory. Set ', ...
             'para.full_racc.solver = ''mosek''.']);
    end
    if strcmpi(opts.formulation, 'raw') && ...
            max(MB * L + 1, MD * L + 1) > opts.raw_max_dimension
        error('FullRACC:RawProblemTooLarge', ...
            ['The raw LMI dimension is %d, exceeding raw_max_dimension=%d. ', ...
             'Use the mathematically equivalent decomposed formulation.'], ...
            max(MB * L + 1, MD * L + 1), opts.raw_max_dimension);
    end
    if ~isempty(opts.checkpoint_dir) && ~exist(opts.checkpoint_dir, 'dir')
        mkdir(opts.checkpoint_dir);
    end

    w_all = zeros(L, num_freqs);
    statuses = repmat({''}, 1, num_freqs);
    relaxation_alpha = nan(1, num_freqs);
    relaxation_ac_db = nan(1, num_freqs);
    recovered_alpha = nan(1, num_freqs);
    recovered_ac_db = nan(1, num_freqs);
    rank_fraction = nan(1, num_freqs);
    rank_ratio = nan(1, num_freqs);
    solve_time_seconds = nan(1, num_freqs);
    etaB_all = nan(1, num_freqs);
    etaD_all = nan(1, num_freqs);
    deltaB_all = nan(1, num_freqs);
    deltaD_all = nan(1, num_freqs);
    bisection_iterations = zeros(1, num_freqs);
    sB_all = nan(1, num_freqs);
    sD_all = nan(1, num_freqs);
    W_all = nan(L, L, num_freqs);
    recovery_source = repmat({''}, 1, num_freqs);
    recovery_iterations = zeros(1, num_freqs);
    normalization_multiplier = nan(1, num_freqs);
    solution_method_used = repmat({''}, 1, num_freqs);
    normalization_attempt_multipliers = cell(1, num_freqs);
    normalization_attempt_statuses = cell(1, num_freqs);

    selected_indices = fre_indices(:).';
    frequency_results = cell(1, numel(selected_indices));
    if opts.parallel_workers > 1 && numel(selected_indices) > 1
        local_ensure_parallel_pool(opts.parallel_workers);
        parfor (position = 1:numel(selected_indices), opts.parallel_workers)
            i = selected_indices(position);
            frequency_results{position} = local_solve_frequency( ...
                HB_ctrl(:, :, i), HD_ctrl(:, :, i), para, opts, i, ...
                num_freqs);
        end
    else
        for position = 1:numel(selected_indices)
            i = selected_indices(position);
            frequency_results{position} = local_solve_frequency( ...
                HB_ctrl(:, :, i), HD_ctrl(:, :, i), para, opts, i, ...
                num_freqs);
        end
    end

    for position = 1:numel(selected_indices)
        i = selected_indices(position);
        item = frequency_results{position};
        w_all(:, i) = item.w;
        statuses{i} = item.status;
        relaxation_alpha(i) = item.relaxation_alpha;
        relaxation_ac_db(i) = item.relaxation_ac_db;
        recovered_alpha(i) = item.recovered_alpha;
        recovered_ac_db(i) = item.recovered_ac_db;
        rank_fraction(i) = item.rank_fraction;
        rank_ratio(i) = item.rank_ratio;
        solve_time_seconds(i) = item.solve_time_seconds;
        etaB_all(i) = item.etaB;
        etaD_all(i) = item.etaD;
        deltaB_all(i) = item.deltaB;
        deltaD_all(i) = item.deltaD;
        bisection_iterations(i) = item.iterations;
        sB_all(i) = item.sB;
        sD_all(i) = item.sD;
        W_all(:, :, i) = item.W;
        recovery_source{i} = item.recovery_source;
        recovery_iterations(i) = item.recovery_iterations;
        if isfield(item, 'normalization_multiplier')
            normalization_multiplier(i) = item.normalization_multiplier;
        end
        if isfield(item, 'solution_method_used')
            solution_method_used{i} = item.solution_method_used;
        end
        if isfield(item, 'attempt_multipliers')
            normalization_attempt_multipliers{i} = item.attempt_multipliers;
        end
        if isfield(item, 'attempt_statuses')
            normalization_attempt_statuses{i} = item.attempt_statuses;
        end
    end

    res = struct();
    res.w = w_all;
    res.scale = 1;
    res.status = statuses;
    res.formulation = opts.formulation;
    res.solver = opts.solver;
    res.eta = struct('B', etaB_all, 'D', etaD_all);
    res.diagonal_loading = struct('B', deltaB_all, 'D', deltaD_all, ...
        'ratio', opts.diagonal_loading_ratio);
    res.relaxation_alpha = relaxation_alpha;
    res.relaxation_ac_db = relaxation_ac_db;
    res.recovered_alpha = recovered_alpha;
    res.recovered_ac_db = recovered_ac_db;
    res.rank_fraction = rank_fraction;
    res.rank_ratio = rank_ratio;
    res.solve_time_seconds = solve_time_seconds;
    res.bisection_iterations = bisection_iterations;
    res.sB = sB_all;
    res.sD = sD_all;
    res.W = W_all;
    res.recovery_source = recovery_source;
    res.recovery_iterations = recovery_iterations;
    res.normalization_multiplier = normalization_multiplier;
    res.solution_method_used = solution_method_used;
    res.normalization_attempt_multipliers = ...
        normalization_attempt_multipliers;
    res.normalization_attempt_statuses = normalization_attempt_statuses;
    res.options = opts;
end

function item = local_solve_frequency(HB, HD, para, opts, i, num_freqs)
    etaB = local_resolve_eta(para, 'B', i, HB);
    etaD = local_resolve_eta(para, 'D', i, HD);
    deltaB = opts.diagonal_loading_ratio * norm(HB, 2)^2;
    deltaD = opts.diagonal_loading_ratio * norm(HD, 2)^2;
    signature = local_checkpoint_signature( ...
        HB, HD, etaB, etaD, deltaB, deltaD, opts, i);
    checkpoint_file = local_checkpoint_file(opts, i);

    cached_item = [];
    if opts.resume_from_checkpoints && ~isempty(checkpoint_file) && ...
            exist(checkpoint_file, 'file')
        saved = load(checkpoint_file, 'checkpoint');
        if isfield(saved, 'checkpoint') && ...
                isequaln(saved.checkpoint.signature, signature)
            cached_item = saved.checkpoint.item;
            cached_item = local_upgrade_cached_item(cached_item);
            if local_status_is_solved(cached_item.status) || ...
                    ~opts.retry_unsolved_checkpoints
                item = cached_item;
                fprintf(['Full-RACC frequency %d/%d restored from ', ...
                    'checkpoint.\n'], i, num_freqs);
                return;
            end
            fprintf(['Full-RACC frequency %d/%d checkpoint status is ', ...
                '%s; retrying with alternative normalization.\n'], ...
                i, num_freqs, cached_item.status);
        end
    end

    base_data_scale = max(norm([HB; HD], 'fro'), 1e-12);
    elapsed = 0;
    solution = [];
    selected_multiplier = NaN;
    method_used = opts.solution_method;
    attempt_multipliers = zeros(1, 0);
    attempt_statuses = cell(1, 0);

    if strcmpi(opts.solution_method, 'normalized')
        multipliers = opts.normalization_multipliers(:).';
        for attempt = 1:numel(multipliers)
            multiplier = multipliers(attempt);
            attempt_multipliers(end + 1) = multiplier; %#ok<AGROW>

            % The previous official run already tried multiplier 1. Reuse
            % that failed status instead of spending another five minutes
            % reproducing the same numerical result.
            if ~isempty(cached_item) && numel(multipliers) > 1 && ...
                    abs(multiplier - 1) <= eps
                attempt_statuses{end + 1} = cached_item.status; %#ok<AGROW>
                continue;
            end

            [HBs, HDs, etaBs, etaDs, deltaBs, deltaDs] = ...
                local_scale_problem(HB, HD, etaB, etaD, deltaB, deltaD, ...
                base_data_scale, multiplier);
            attempt_opts = opts;
            if ~isempty(opts.normalization_sD)
                attempt_opts.normalization_sD = ...
                    opts.normalization_sD / multiplier^2;
            end
            attempt_timer = tic;
            candidate = local_optimize_normalized( ...
                HBs, HDs, etaBs, etaDs, deltaBs, deltaDs, attempt_opts);
            elapsed = elapsed + toc(attempt_timer);
            attempt_statuses{end + 1} = candidate.status; %#ok<AGROW>
            fprintf(['Full-RACC frequency %d/%d normalization multiplier ', ...
                '%g returned %s.\n'], i, num_freqs, multiplier, ...
                candidate.status);
            solution = candidate;
            if local_status_is_solved(candidate.status)
                selected_multiplier = multiplier;
                break;
            end
        end

        if (isempty(solution) || ~local_status_is_solved(solution.status)) && ...
                opts.fallback_to_bisection
            [HBs, HDs, etaBs, etaDs, deltaBs, deltaDs] = ...
                local_scale_problem(HB, HD, etaB, etaD, deltaB, deltaD, ...
                base_data_scale, 1);
            fallback_timer = tic;
            solution = local_bisect( ...
                HBs, HDs, etaBs, etaDs, deltaBs, deltaDs, opts);
            elapsed = elapsed + toc(fallback_timer);
            method_used = 'bisection fallback';
            fprintf(['Full-RACC frequency %d/%d normalization retries ', ...
                'failed; bisection fallback returned %s.\n'], ...
                i, num_freqs, solution.status);
        end
    elseif strcmpi(opts.solution_method, 'bisection')
        [HBs, HDs, etaBs, etaDs, deltaBs, deltaDs] = ...
            local_scale_problem(HB, HD, etaB, etaD, deltaB, deltaD, ...
            base_data_scale, 1);
        freq_timer = tic;
        solution = local_bisect( ...
            HBs, HDs, etaBs, etaDs, deltaBs, deltaDs, opts);
        elapsed = toc(freq_timer);
    else
        error('FullRACC:UnknownSolutionMethod', ...
            'Unknown solution method "%s".', opts.solution_method);
    end

    previous_elapsed = 0;
    if ~isempty(cached_item) && isfield(cached_item, 'solve_time_seconds')
        previous_elapsed = cached_item.solve_time_seconds;
    end

    item = struct('w', zeros(size(HB, 2), 1), 'status', solution.status, ...
        'relaxation_alpha', NaN, 'relaxation_ac_db', NaN, ...
        'recovered_alpha', NaN, 'recovered_ac_db', NaN, ...
        'rank_fraction', NaN, 'rank_ratio', NaN, ...
        'solve_time_seconds', previous_elapsed + elapsed, ...
        'latest_retry_time_seconds', elapsed, ...
        'etaB', etaB, 'etaD', etaD, ...
        'deltaB', deltaB, 'deltaD', deltaD, ...
        'iterations', solution.iterations, 'sB', solution.sB, ...
        'sD', solution.sD, 'W', nan(size(HB, 2)), ...
        'recovery_source', '', 'recovery_iterations', 0, ...
        'normalization_multiplier', selected_multiplier, ...
        'solution_method_used', method_used, ...
        'attempt_multipliers', attempt_multipliers, ...
        'attempt_statuses', {attempt_statuses});

    if isempty(solution.W)
        warning('FullRACC:FrequencyUnsolved', ...
            ['Full-RACC did not return a feasible solution at frequency ', ...
             'index %d (%s).'], i, solution.status);
    else
        [w, recovery] = full_racc_recover_filter( ...
            solution.W, HB, HD, etaB, etaD, opts, i, deltaB, deltaD);
        [ac_db_i, alpha_i] = full_racc_worst_case_contrast( ...
            w, HB, HD, etaB, etaD, deltaB, deltaD);
        item.w = w;
        item.relaxation_alpha = solution.alpha;
        item.relaxation_ac_db = solution.alpha_db;
        item.recovered_alpha = alpha_i;
        item.recovered_ac_db = ac_db_i;
        item.rank_fraction = recovery.rank_fraction;
        item.rank_ratio = recovery.rank_ratio;
        item.W = solution.W;
        item.recovery_source = recovery.source;
        item.recovery_iterations = recovery.refinement_iterations;
        fprintf(['Full-RACC frequency %d/%d: relaxation %.2f dB, ', ...
            'recovered %.2f dB, rank fraction %.4f, %.1f s.\n'], ...
            i, num_freqs, solution.alpha_db, ac_db_i, ...
            recovery.rank_fraction, item.solve_time_seconds);
    end

    if ~isempty(checkpoint_file)
        checkpoint = struct('signature', signature, 'item', item); %#ok<NASGU>
        save(checkpoint_file, 'checkpoint', '-v7.3');
    end
end

function [HBs, HDs, etaBs, etaDs, deltaBs, deltaDs] = ...
        local_scale_problem(HB, HD, etaB, etaD, deltaB, deltaD, ...
        base_data_scale, multiplier)
% RACC_PM_Sub uses the same numerical strategy: scale H and its uncertainty
% radii together, retrying several multipliers. The physical robust problem
% is unchanged because every zone-energy quantity scales consistently.
    data_scale = base_data_scale * multiplier;
    HBs = HB / data_scale;
    HDs = HD / data_scale;
    etaBs = etaB / data_scale;
    etaDs = etaD / data_scale;
    deltaBs = deltaB / data_scale^2;
    deltaDs = deltaD / data_scale^2;
end

function item = local_upgrade_cached_item(item)
% Checkpoints written before normalization-retry diagnostics were added all
% used the baseline multiplier 1 and the normalized one-SDP method.
    if ~isfield(item, 'normalization_multiplier')
        item.normalization_multiplier = 1;
    end
    if ~isfield(item, 'solution_method_used')
        item.solution_method_used = 'normalized';
    end
    if ~isfield(item, 'attempt_multipliers')
        item.attempt_multipliers = 1;
    end
    if ~isfield(item, 'attempt_statuses')
        item.attempt_statuses = {item.status};
    end
    if ~isfield(item, 'latest_retry_time_seconds')
        item.latest_retry_time_seconds = 0;
    end
end

function solution = local_optimize_normalized( ...
        HB, HD, etaB, etaD, deltaBLoading, deltaDLoading, opts)
% Replace trace(W)=1 by a fixed positive s_D and maximize s_B. All robust
% LMIs are positively homogeneous in (W,s_B,s_D,lambda_B,lambda_D), so any
% positive fixed value is exact. A data-scaled s_D keeps the variables near
% their natural magnitudes when eta_D is small, while s_D=1 can make an
% ill-conditioned problem appear spuriously unbounded to a conic solver.
    MB = size(HB, 1);
    MD = size(HD, 1);
    L = size(HB, 2);
    sDReference = local_normalization_sD(HD, etaD, deltaDLoading, opts);
    traceUpper = Inf;
    if etaD > 0
        traceUpper = min(traceUpper, L * sDReference / etaD^2);
    end
    if deltaDLoading > 0
        traceUpper = min(traceUpper, sDReference / deltaDLoading);
    end
    brightEnergyCoefficient = norm(HB, 2)^2 + deltaBLoading;

    if strcmpi(opts.formulation, 'decomposed')
        cvx_begin sdp quiet
            cvx_solver(opts.solver)
            if ~isempty(opts.mosek_threads)
                cvx_solver_settings('MSK_IPAR_NUM_THREADS', opts.mosek_threads)
            end
            variable W(L, L) hermitian semidefinite
            variable sB nonnegative
            variable lambdaB nonnegative
            variable lambdaD nonnegative
            variable deltaB(MB) nonnegative
            variable deltaD(MD) nonnegative
            expression UB(MB, L)
            expression UD(MD, L)
            expression cB
            expression cD
            UB = HB * W;
            UD = HD * W;
            cB = real(trace(HB * W * HB')) + ...
                deltaBLoading * real(trace(W));
            cD = real(trace(HD * W * HD')) + ...
                deltaDLoading * real(trace(W));
            maximize(sB)
            subject to
                if etaD > 0
                    % These bounds are already implied by the dark-zone
                    % LMI and fixed s_D, but stating them explicitly prevents
                    % large ill-conditioned problems from being reported
                    % as spuriously unbounded by the conic solver.
                    lambdaD <= sDReference / etaD^2;
                end
                if isfinite(traceUpper)
                    real(trace(W)) <= traceUpper;
                    sB <= brightEnergyCoefficient * traceUpper;
                end
                for m = 1:MB
                    uBm = UB(m, :).';
                    [W.' + lambdaB * eye(L), uBm; ...
                     uBm', deltaB(m)] >= 0;
                end
                sum(deltaB) <= cB - sB - lambdaB * etaB^2;
                for m = 1:MD
                    uDm = UD(m, :).';
                    [lambdaD * eye(L) - W.', -uDm; ...
                     -uDm', deltaD(m)] >= 0;
                end
                sum(deltaD) <= sDReference - cD - lambdaD * etaD^2;
        cvx_end
    elseif strcmpi(opts.formulation, 'raw')
        cvx_begin sdp quiet
            cvx_solver(opts.solver)
            if startsWith(lower(opts.solver), 'mosek') && ~isempty(opts.mosek_threads)
                cvx_solver_settings('MSK_IPAR_NUM_THREADS', opts.mosek_threads)
            end
            variable W(L, L) hermitian semidefinite
            variable sB nonnegative
            variable lambdaB nonnegative
            variable lambdaD nonnegative
            expression QB(MB * L, MB * L)
            expression QD(MD * L, MD * L)
            expression uB(MB * L)
            expression uD(MD * L)
            expression cB
            expression cD
            QB = kron(W.', eye(MB));
            QD = kron(W.', eye(MD));
            uB = reshape(HB * W, MB * L, 1);
            uD = reshape(HD * W, MD * L, 1);
            cB = real(trace(HB * W * HB')) + ...
                deltaBLoading * real(trace(W));
            cD = real(trace(HD * W * HD')) + ...
                deltaDLoading * real(trace(W));
            maximize(sB)
            subject to
                if etaD > 0
                    lambdaD <= sDReference / etaD^2;
                end
                if isfinite(traceUpper)
                    real(trace(W)) <= traceUpper;
                    sB <= brightEnergyCoefficient * traceUpper;
                end
                [QB + lambdaB * eye(MB * L), uB; ...
                 uB', cB - sB - lambdaB * etaB^2] >= 0;
                [-QD + lambdaD * eye(MD * L), -uD; ...
                 -uD', sDReference - cD - lambdaD * etaD^2] >= 0;
        cvx_end
    else
        error('FullRACC:UnknownFormulation', ...
            'Unknown formulation "%s".', opts.formulation);
    end

    feasible = strcmp(cvx_status, 'Solved') || ...
        strcmp(cvx_status, 'Inaccurate/Solved');
    if feasible
        alpha = (MD / MB) * full(sB) / sDReference;
        if alpha > 0
            alpha_db = 10 * log10(alpha);
        else
            alpha_db = -Inf;
        end
        solution = struct('W', full((W + W') / 2), ...
            'alpha', alpha, 'alpha_db', alpha_db, 'status', cvx_status, ...
            'iterations', 1, 'sB', full(sB), 'sD', sDReference);
    else
        solution = struct('W', [], 'alpha', NaN, 'alpha_db', NaN, ...
            'status', cvx_status, 'iterations', 1, 'sB', NaN, ...
            'sD', sDReference);
    end
end

function solution = local_bisect( ...
        HB, HD, etaB, etaD, deltaBLoading, deltaDLoading, opts)
% For fixed alpha, local_feasibility is an SDP. The returned value is
%
%   alpha_SDR^* = sup {alpha >= 0 : F_SDR(alpha) is feasible}.
%
% The lower/upper search is performed in dB only for numerical convenience.
    [MB, L] = size(HB);
    MD = size(HD, 1);
    RB = (HB' * HB + (HB' * HB)') / 2 + ...
        deltaBLoading * eye(L);
    RD = (HD' * HD + (HD' * HD)') / 2 + ...
        deltaDLoading * eye(L);

    rd_floor = max(1e-10 * max(real(trace(RD)) / L, 1), eps);
    [~, D] = eig(RB, RD + rd_floor * eye(L));
    nominal_alpha = (MD / MB) * max(real(diag(D)));
    nominal_db = 10 * log10(max(nominal_alpha, realmin));

    lower_db = opts.alpha_floor_db;
    if isempty(opts.alpha_upper_db)
        upper_db = min(nominal_db + opts.upper_margin_db, opts.maximum_upper_db);
        upper_db = max(upper_db, lower_db + 1);
    else
        upper_db = opts.alpha_upper_db;
    end

    [upper_feasible, upper_state] = local_feasibility( ...
        HB, HD, etaB, etaD, deltaBLoading, deltaDLoading, ...
        10^(upper_db / 10), opts);
    while upper_feasible && upper_db < opts.maximum_upper_db
        lower_db = upper_db;
        lower_state = upper_state;
        upper_db = min(upper_db + 20, opts.maximum_upper_db);
        [upper_feasible, upper_state] = local_feasibility( ...
            HB, HD, etaB, etaD, deltaBLoading, deltaDLoading, ...
            10^(upper_db / 10), opts);
    end

    if ~exist('lower_state', 'var')
        [lower_feasible, lower_state] = local_feasibility( ...
            HB, HD, etaB, etaD, deltaBLoading, deltaDLoading, ...
            10^(lower_db / 10), opts);
        if ~lower_feasible
            solution = struct('W', [], 'alpha', NaN, 'alpha_db', NaN, ...
                'status', ['Lower-bound infeasible: ' lower_state.status], ...
                'iterations', 0, 'sB', NaN, 'sD', NaN);
            return;
        end
    end

    if upper_feasible
        solution = local_pack_solution(upper_state, upper_db, ...
            'Upper search limit remained feasible', 0);
        return;
    end

    iteration = 0;
    while (upper_db - lower_db) > opts.bisection_tolerance_db && ...
            iteration < opts.max_bisection_iterations
        iteration = iteration + 1;
        candidate_db = (lower_db + upper_db) / 2;
        [is_feasible, state] = local_feasibility( ...
            HB, HD, etaB, etaD, deltaBLoading, deltaDLoading, ...
            10^(candidate_db / 10), opts);
        if is_feasible
            lower_db = candidate_db;
            lower_state = state;
        else
            upper_db = candidate_db;
        end
    end

    solution = local_pack_solution(lower_state, lower_db, ...
        lower_state.status, iteration);
end

function solution = local_pack_solution(state, alpha_db, status, iterations)
    solution = struct();
    solution.W = state.W;
    solution.alpha_db = alpha_db;
    solution.alpha = 10^(alpha_db / 10);
    solution.status = status;
    solution.iterations = iterations;
    solution.sB = state.sB;
    solution.sD = state.sD;
end

function [feasible, state] = local_feasibility( ...
        HB, HD, etaB, etaD, deltaBLoading, deltaDLoading, alpha, opts)
    MB = size(HB, 1);
    MD = size(HD, 1);
    L = size(HB, 2);

    if strcmpi(opts.formulation, 'decomposed')
        cvx_begin sdp quiet
            cvx_solver(opts.solver)
            if ~isempty(opts.mosek_threads)
                cvx_solver_settings('MSK_IPAR_NUM_THREADS', opts.mosek_threads)
            end
            variable W(L, L) hermitian semidefinite
            variable sB nonnegative
            variable sD nonnegative
            variable lambdaB nonnegative
            variable lambdaD nonnegative
            variable deltaB(MB) nonnegative
            variable deltaD(MD) nonnegative
            expression UB(MB, L)
            expression UD(MD, L)
            expression cB
            expression cD
            UB = HB * W;
            UD = HD * W;
            cB = real(trace(HB * W * HB')) + ...
                deltaBLoading * real(trace(W));
            cD = real(trace(HD * W * HD')) + ...
                deltaDLoading * real(trace(W));
            minimize(0)
            subject to
                real(trace(W)) == 1;
                MD * sB - alpha * MB * sD >= 0;
                for m = 1:MB
                    % The ordinary transpose is required by
                    % x=vec(Delta H), Q=W^T kron I_M.
                    uBm = UB(m, :).';
                    [W.' + lambdaB * eye(L), uBm; ...
                     uBm', deltaB(m)] >= 0;
                end
                sum(deltaB) <= cB - sB - lambdaB * etaB^2;
                for m = 1:MD
                    uDm = UD(m, :).';
                    [lambdaD * eye(L) - W.', -uDm; ...
                     -uDm', deltaD(m)] >= 0;
                end
                sum(deltaD) <= sD - cD - lambdaD * etaD^2;
        cvx_end
    elseif strcmpi(opts.formulation, 'raw')
        cvx_begin sdp quiet
            cvx_solver(opts.solver)
            if startsWith(lower(opts.solver), 'mosek') && ~isempty(opts.mosek_threads)
                cvx_solver_settings('MSK_IPAR_NUM_THREADS', opts.mosek_threads)
            end
            variable W(L, L) hermitian semidefinite
            variable sB nonnegative
            variable sD nonnegative
            variable lambdaB nonnegative
            variable lambdaD nonnegative
            expression QB(MB * L, MB * L)
            expression QD(MD * L, MD * L)
            expression uB(MB * L)
            expression uD(MD * L)
            expression cB
            expression cD
            QB = kron(W.', eye(MB));
            QD = kron(W.', eye(MD));
            uB = reshape(HB * W, MB * L, 1);
            uD = reshape(HD * W, MD * L, 1);
            cB = real(trace(HB * W * HB')) + ...
                deltaBLoading * real(trace(W));
            cD = real(trace(HD * W * HD')) + ...
                deltaDLoading * real(trace(W));
            minimize(0)
            subject to
                real(trace(W)) == 1;
                MD * sB - alpha * MB * sD >= 0;
                [QB + lambdaB * eye(MB * L), uB; ...
                 uB', cB - sB - lambdaB * etaB^2] >= 0;
                [-QD + lambdaD * eye(MD * L), -uD; ...
                 -uD', sD - cD - lambdaD * etaD^2] >= 0;
        cvx_end
    else
        error('FullRACC:UnknownFormulation', ...
            'Unknown formulation "%s".', opts.formulation);
    end

    feasible = strcmp(cvx_status, 'Solved') || ...
        strcmp(cvx_status, 'Inaccurate/Solved');
    state = struct('status', cvx_status, 'W', [], 'sB', NaN, 'sD', NaN);
    if feasible
        state.W = full((W + W') / 2);
        state.sB = full(sB);
        state.sD = full(sD);
    end
end

function opts = local_options(para)
    opts = struct( ...
        'solution_method', 'normalized', ...
        'formulation', 'decomposed', ...
        'solver', 'mosek', ...
        'mosek_threads', 4, ...
        'bisection_tolerance_db', 0.10, ...
        'max_bisection_iterations', 14, ...
        'alpha_floor_db', -80, ...
        'alpha_upper_db', [], ...
        'upper_margin_db', 10, ...
        'maximum_upper_db', 160, ...
        'randomization_trials', 300, ...
        'random_seed', 20260820, ...
        'recovery_method', 'paper', ...
        'refinement_iterations', 300, ...
        'refinement_seeds', 8, ...
        'refinement_tolerance', 1e-8, ...
        'parallel_workers', 1, ...
        'checkpoint_dir', '', ...
        'resume_from_checkpoints', true, ...
        'retry_unsolved_checkpoints', true, ...
        'normalization_multipliers', 1, ...
        'fallback_to_bisection', false, ...
        'normalization_sD', [], ...
        'diagonal_loading_ratio', 0, ...
        'raw_max_dimension', 800);
    if isfield(para, 'full_racc')
        names = fieldnames(opts);
        for k = 1:numel(names)
            if isfield(para.full_racc, names{k})
                opts.(names{k}) = para.full_racc.(names{k});
            end
        end
    end
    opts.formulation = char(string(opts.formulation));
    opts.solution_method = char(string(opts.solution_method));
    opts.solver = char(string(opts.solver));
    opts.recovery_method = char(string(opts.recovery_method));
    opts.checkpoint_dir = char(string(opts.checkpoint_dir));
    validateattributes(opts.diagonal_loading_ratio, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(opts.parallel_workers, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(opts.randomization_trials, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative'});
    validateattributes(opts.resume_from_checkpoints, ...
        {'logical', 'numeric'}, {'scalar'});
    opts.resume_from_checkpoints = logical(opts.resume_from_checkpoints);
    validateattributes(opts.retry_unsolved_checkpoints, ...
        {'logical', 'numeric'}, {'scalar'});
    opts.retry_unsolved_checkpoints = ...
        logical(opts.retry_unsolved_checkpoints);
    validateattributes(opts.fallback_to_bisection, ...
        {'logical', 'numeric'}, {'scalar'});
    opts.fallback_to_bisection = logical(opts.fallback_to_bisection);
    validateattributes(opts.normalization_multipliers, {'numeric'}, ...
        {'vector', 'real', 'finite', 'positive'});
    opts.normalization_multipliers = ...
        unique(opts.normalization_multipliers(:).', 'stable');
end

function signature = local_checkpoint_signature( ...
        HB, HD, etaB, etaD, deltaB, deltaD, opts, frequency_index)
    signature = struct( ...
        'version', 1, ...
        'frequency_index', frequency_index, ...
        'size_HB', size(HB), 'size_HD', size(HD), ...
        'HB_power', norm(HB, 'fro')^2, ...
        'HD_power', norm(HD, 'fro')^2, ...
        'HB_sum', sum(HB(:)), 'HD_sum', sum(HD(:)), ...
        'etaB', etaB, 'etaD', etaD, ...
        'deltaB', deltaB, 'deltaD', deltaD, ...
        'solution_method', opts.solution_method, ...
        'formulation', opts.formulation, 'solver', opts.solver, ...
        'recovery_method', opts.recovery_method, ...
        'randomization_trials', opts.randomization_trials, ...
        'random_seed', opts.random_seed);
end

function checkpoint_file = local_checkpoint_file(opts, frequency_index)
    if isempty(opts.checkpoint_dir)
        checkpoint_file = '';
    else
        checkpoint_file = fullfile(opts.checkpoint_dir, sprintf( ...
            'full_wcracc_frequency_%04d.mat', frequency_index));
    end
end

function local_ensure_parallel_pool(requested_workers)
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', requested_workers);
    elseif pool.NumWorkers < requested_workers
        warning('FullRACC:ExistingPoolSmallerThanRequested', ...
            ['The existing pool has %d workers; requested %d. The existing ', ...
             'pool is retained to avoid discarding active worker state.'], ...
            pool.NumWorkers, requested_workers);
    end
end

function solved = local_status_is_solved(status)
    solved = strcmp(status, 'Solved') || strcmp(status, 'Inaccurate/Solved');
end

function sDReference = local_normalization_sD( ...
        HD, etaD, deltaDLoading, opts)
% Use the nominal mean dark-zone energy generated by an isotropic unit-trace
% matrix as a dimensionally consistent normalization. eta_D covers the
% nearly-zero nominal dark-zone edge case. Any positive reference is exact
% because the normalized robust contrast SDP is positively homogeneous.
    if isempty(opts.normalization_sD)
        L = size(HD, 2);
        nominal_mean_energy = norm(HD, 'fro')^2 / L + deltaDLoading;
        sDReference = max([nominal_mean_energy, etaD^2, eps]);
    else
        sDReference = opts.normalization_sD;
    end
    if ~isnumeric(sDReference) || ~isscalar(sDReference) || ...
            ~isreal(sDReference) || ~isfinite(sDReference) || ...
            sDReference <= 0
        error('FullRACC:InvalidNormalizationSD', ...
            'normalization_sD must be a finite positive scalar.');
    end
end

function eta = local_resolve_eta(para, zone, freq_idx, H)
    if isfield(para, 'full_racc') && isfield(para.full_racc, 'eta') && ...
            isfield(para.full_racc.eta, zone)
        eta_source = para.full_racc.eta.(zone);
    elseif isfield(para, 'epsilon') && isfield(para.epsilon, zone)
        eta_source = para.epsilon.(zone);
    elseif isfield(para, 'full_racc') && isfield(para.full_racc, 'relative_eta')
        relative_source = para.full_racc.relative_eta;
        if isstruct(relative_source)
            relative_source = relative_source.(zone);
        end
        eta = local_pick_value(relative_source, freq_idx) * norm(H, 'fro');
        local_validate_eta(eta, zone, freq_idx);
        return;
    else
        error('FullRACC:MissingUncertaintyRadius', ...
            ['Provide para.full_racc.eta.%s, para.epsilon.%s, or ', ...
             'para.full_racc.relative_eta.'], zone, zone);
    end
    eta = local_pick_value(eta_source, freq_idx);
    local_validate_eta(eta, zone, freq_idx);
end

function value = local_pick_value(source, freq_idx)
    if isscalar(source)
        value = source;
    elseif numel(source) >= freq_idx
        value = source(freq_idx);
    else
        error('FullRACC:UncertaintyLengthMismatch', ...
            'The uncertainty-radius vector does not contain frequency index %d.', freq_idx);
    end
end

function local_validate_eta(eta, zone, freq_idx)
    if ~isnumeric(eta) || ~isscalar(eta) || ~isreal(eta) || ...
            ~isfinite(eta) || eta < 0
        error('FullRACC:InvalidUncertaintyRadius', ...
            'eta_%s at frequency index %d must be a finite nonnegative scalar.', ...
            zone, freq_idx);
    end
end
