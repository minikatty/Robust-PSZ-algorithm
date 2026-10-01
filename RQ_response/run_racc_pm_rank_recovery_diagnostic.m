function results = run_racc_pm_rank_recovery_diagnostic(user_config)
%RUN_RACC_PM_RANK_RECOVERY_DIAGNOSTIC Quantify RACC-PM SDR recovery.
%
% This reviewer-facing diagnostic reproduces the maintained no-suffix
% RACC_PM_Sub LMI, but retains W_tilde so that the relaxation rank and the
% direct relaxation-to-recovered-filter AC/NSRE gaps can be reported. The
% default audit uses the final rho=10, nu=0.01 setting over 100--4000 Hz.
% It is deliberately separate from the production solver and must not be run
% while another large MOSEK job is active.

    if nargin < 1 || isempty(user_config)
        user_config = struct();
    end
    response_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(response_dir);
    addpath(fullfile(project_root, 'src', 'algorithm'));
    addpath(fullfile(project_root, 'src', 'utils'));
    addpath(fullfile(project_root, 'src', 'evaluations'), '-begin');
    config = local_default_config(response_dir);
    config = local_apply_overrides(config, user_config);
    local_validate_config(config);
    rng(config.random_seed, 'twister');

    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    if strlength(string(config.resume_output_dir)) > 0
        output_dir = char(config.resume_output_dir);
    else
        output_dir = fullfile(config.output_root, ['run_' timestamp]);
    end
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    data_dir = fullfile(project_root, 'data', 'SimulateRIR', 'temperature');
    metadata = load(fullfile(data_dir, 'para.mat'), 'para');
    all_frequencies_hz = metadata.para.freq_params.target_freqs(:).';
    frequency_indices = local_frequency_indices(all_frequencies_hz, ...
        config.requested_frequencies_hz, config.frequency_range_hz);
    frequencies_hz = all_frequencies_hz(frequency_indices);
    nominal_file = get_data_filename(data_dir, 'temperature', ...
        config.nominal_temperature_c);
    nominal = load(nominal_file, 'ATF_BZ', 'ATF_DZ');
    HB = nominal.ATF_BZ.ctrl(:, :, frequency_indices);
    HD = nominal.ATF_DZ.ctrl(:, :, frequency_indices);
    desired_data = load(fullfile(project_root, 'data', ...
        'ATF_desired_plane.mat'), 'ATF_desired_plane');
    desired = desired_data.ATF_desired_plane(:, frequency_indices);
    [MB, L, num_frequencies] = size(HB);
    MD = size(HD, 1);

    [epsilon_B, epsilon_D] = local_radii(HB, HD, ...
        config.relative_uncertainty_level);
    [reference_B, reference_D] = local_radii(HB, HD, ...
        config.alpha_reference_relative_level);
    reference_para.full_racc.eta = struct('B', reference_B, ...
        'D', reference_D);
    reference_para.full_racc.diagonal_loading_ratio = 0;
    reference = NoCT_WCRACC(HB, HD, reference_para);
    reference_ac_db = nan(1, num_frequencies);
    for k = 1:num_frequencies
        reference_ac_db(k) = local_ac_db(reference.w(:, k), ...
            HB(:, :, k), HD(:, :, k));
    end
    alpha = 10.^((reference_ac_db - config.alpha_margin_db) / 10);
    gamma = config.mu + config.rho * alpha;

    rows = cell(num_frequencies, 1);
    solver_outputs = cell(num_frequencies, 1);
    for k = 1:num_frequencies
        checkpoint_file = fullfile(output_dir, sprintf( ...
            'frequency_%04dHz.mat', round(frequencies_hz(k))));
        if exist(checkpoint_file, 'file')
            checkpoint = load(checkpoint_file, 'output', 'recovery', 'row');
            if isfield(checkpoint, 'row') && isfield(checkpoint, 'output')
                rows{k} = checkpoint.row;
                solver_outputs{k} = checkpoint.output;
                fprintf('Rank diagnostic %.0f Hz (%d/%d) restored.\n', ...
                    frequencies_hz(k), k, num_frequencies);
                continue;
            end
        end
        fprintf('Rank diagnostic %.0f Hz (%d/%d).\n', ...
            frequencies_hz(k), k, num_frequencies);
        H_energy = max(norm([HB(:, :, k); HD(:, :, k)], 'fro'), ...
            realmin);
        base_scale = H_energy * sqrt(config.rho);
        output = [];
        attempt_status = strings(numel(config.normalization_multipliers), 1);
        for attempt = 1:numel(config.normalization_multipliers)
            multiplier = config.normalization_multipliers(attempt);
            scale = base_scale * multiplier;
            output = local_solve_lmi(HB(:, :, k) / scale, ...
                HD(:, :, k) / scale, desired(:, k) / scale, ...
                epsilon_B(k) / scale, epsilon_D(k) / scale, ...
                config.rho, gamma(k), config.power_bound, MB, MD, L);
            attempt_status(attempt) = string(output.status);
            if local_status_is_solved(output.status)
                output.normalization_multiplier = multiplier;
                break;
            end
        end
        output.attempt_status = attempt_status;
        if ~local_status_is_solved(output.status)
            warning('RACCPMRankDiagnostic:Unsolved', ...
                'No solved LMI at %.0f Hz.', frequencies_hz(k));
            row = local_unsolved_row(frequencies_hz(k), ...
                epsilon_B(k), epsilon_D(k), alpha(k), output);
            rows{k} = row;
            solver_outputs{k} = output;
            recovery = struct();
            save(checkpoint_file, 'output', 'recovery', 'row', '-v7.3');
            continue;
        end

        recovery = local_recover_and_score(output.W_tilde, ...
            HB(:, :, k), HD(:, :, k), desired(:, k), config, ...
            output.scaled_problem);
        output.recovery = recovery;
        solver_outputs{k} = output;
        row = table(frequencies_hz(k), ...
            config.relative_uncertainty_level, epsilon_B(k), epsilon_D(k), ...
            alpha(k), string(output.status), output.solve_time_seconds, ...
            output.normalization_multiplier, recovery.rank_fraction, ...
            recovery.tail_fraction, recovery.rank_ratio, ...
            recovery.second_to_first_ratio, ...
            recovery.relative_frobenius_rank1_error, ...
            recovery.numerical_rank, ...
            string(recovery.production_mode), ...
            recovery.relaxed_ac_db, recovery.relaxed_nsre_db, ...
            recovery.production_ac_db, recovery.production_nsre_db, ...
            recovery.production_minus_relaxed_ac_db, ...
            recovery.production_minus_relaxed_nsre_db, ...
            recovery.principal_ac_db, recovery.principal_nsre_db, ...
            recovery.gaussian_ac_db, recovery.gaussian_nsre_db, ...
            recovery.gaussian_minus_principal_ac_db, ...
            recovery.gaussian_minus_principal_nsre_db, ...
            recovery.random_ac_p10_db, recovery.random_ac_p90_db, ...
            recovery.random_nsre_p10_db, recovery.random_nsre_p90_db, ...
            'VariableNames', local_variable_names());
        rows{k} = row;
        save(checkpoint_file, 'output', 'recovery', 'row', '-v7.3');
    end

    summary = vertcat(rows{:});
    writetable(summary, fullfile(output_dir, ...
        'rank_recovery_summary.csv'));
    band_summary = local_band_summary(summary, config);
    writetable(band_summary, fullfile(output_dir, ...
        'rank_recovery_band_summary.csv'));
    results = struct('config', config, 'output_dir', output_dir, ...
        'frequencies_hz', frequencies_hz, ...
        'frequency_indices', frequency_indices, 'epsilon_B', epsilon_B, ...
        'epsilon_D', epsilon_D, 'alpha', alpha, ...
        'reference_ac_db', reference_ac_db, ...
        'solver_outputs', {solver_outputs}, 'summary', summary, ...
        'band_summary', band_summary);
    save(fullfile(output_dir, 'rank_recovery_results.mat'), ...
        'results', '-v7.3');
    disp(summary);
    disp(band_summary);
    fprintf('Completed rank diagnostic: %s\n', output_dir);
end

function config = local_default_config(response_dir)
    config = struct();
    config.nominal_temperature_c = 22.5;
    config.requested_frequencies_hz = [];
    config.frequency_range_hz = [100, 4000];
    config.relative_uncertainty_level = 0.01;
    config.alpha_reference_relative_level = 0.01;
    config.alpha_margin_db = 0;
    config.rho = 10;
    config.mu = 1;
    config.power_bound = 1e4;
    config.randomization_trials = 1000;
    config.randomize_rank_one_solutions = false;
    config.normalization_multipliers = [0.1, 1, 5, 10, 0.01];
    config.rank_fraction_threshold = 0.98;
    config.rank_tolerance = 1e-7;
    config.random_seed = 20260826;
    config.resume_output_dir = "";
    config.output_root = fullfile(response_dir, ...
        'racc_pm_rank_recovery_results');
end

function config = local_apply_overrides(config, overrides)
    if ~isstruct(overrides) || ~isscalar(overrides)
        error('RACCPMRankDiagnostic:InvalidConfig', ...
            'user_config must be a scalar structure.');
    end
    fields = fieldnames(overrides);
    for k = 1:numel(fields)
        if ~isfield(config, fields{k})
            error('RACCPMRankDiagnostic:UnknownConfigField', ...
                'Unknown configuration field: %s.', fields{k});
        end
        config.(fields{k}) = overrides.(fields{k});
    end
end

function local_validate_config(config)
    if ~isempty(config.requested_frequencies_hz)
        validateattributes(config.requested_frequencies_hz, {'numeric'}, ...
            {'vector', 'real', 'finite', 'positive'});
    end
    validateattributes(config.frequency_range_hz, {'numeric'}, ...
        {'vector', 'numel', 2, 'real', 'finite', 'positive', 'increasing'});
    validateattributes(config.relative_uncertainty_level, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(config.rho, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(config.randomization_trials, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(config.rank_fraction_threshold, {'numeric'}, ...
        {'scalar', 'real', 'finite', '>', 0, '<=', 1});
end

function indices = local_frequency_indices(frequencies_hz, requested_hz, ...
        frequency_range_hz)
    if isempty(requested_hz)
        indices = find(frequencies_hz >= frequency_range_hz(1) & ...
            frequencies_hz <= frequency_range_hz(2));
        if isempty(indices)
            error('RACCPMRankDiagnostic:EmptyFrequencyRange', ...
                'No available frequencies lie in the requested range.');
        end
        return;
    end
    indices = zeros(size(requested_hz));
    for k = 1:numel(requested_hz)
        [difference, indices(k)] = min(abs(frequencies_hz - requested_hz(k)));
        if difference > 1e-9
            error('RACCPMRankDiagnostic:FrequencyUnavailable', ...
                'Requested frequency %.6g Hz is unavailable.', requested_hz(k));
        end
    end
end

function [epsilon_B, epsilon_D] = local_radii(HB, HD, relative_level)
    num_frequencies = size(HB, 3);
    epsilon_B = zeros(1, num_frequencies);
    epsilon_D = zeros(1, num_frequencies);
    for k = 1:num_frequencies
        epsilon_B(k) = relative_level * norm(HB(:, :, k), 'fro');
        epsilon_D(k) = relative_level * norm(HD(:, :, k), 'fro');
    end
end

function output = local_solve_lmi(HB, HD, desired, epsilon_B, epsilon_D, ...
        rho, gamma, power_bound, MB, MD, L)
    sigma_sq = 1e-8;
    timer = tic;
    cvx_clear;
    cvx_begin sdp quiet
        cvx_solver mosek
        variable W_tilde(L+1, L+1) hermitian
        variable tB
        variable tD
        variable tau_B nonnegative
        variable tau_D nonnegative
        variable delta_B(MB)
        variable delta_D(MD)
        minimize(tB + tD)
        subject to
            W_tilde >= 0;
            W_tilde(L+1, L+1) == 1;
            W = W_tilde(1:L, 1:L);
            w = W_tilde(1:L, L+1);
            real(trace(W)) <= power_bound;
            cB = (1-rho) * real(trace(HB * W * HB')) ...
                - 2 * real(desired' * HB * w) + desired' * desired ...
                + sigma_sq - tB;
            sum(delta_B) <= -cB - tau_B * epsilon_B^2;
            for m = 1:MB
                h = HB(m, :);
                p = desired(m);
                u = ((1-rho) * h * W - p * w').';
                [tau_B * eye(L) - (1-rho) * W.', u; ...
                    u', delta_B(m)] >= 0;
            end
            cD = gamma * real(trace(HD * W * HD')) - tD;
            sum(delta_D) <= -cD - tau_D * epsilon_D^2;
            for m = 1:MD
                h = HD(m, :);
                u = (gamma * h * W).';
                [tau_D * eye(L) - gamma * W.', u; ...
                    u', delta_D(m)] >= 0;
            end
    cvx_end
    output = struct();
    output.status = cvx_status;
    output.solve_time_seconds = toc(timer);
    output.objective = cvx_optval;
    output.W_tilde = W_tilde;
    output.scaled_problem = struct('HB', HB, 'HD', HD, ...
        'desired', desired, 'rho', rho, 'gamma', gamma, ...
        'power_bound', power_bound);
end

function recovery = local_recover_and_score(W_tilde, HB, HD, desired, ...
        config, scaled_problem)
    W_tilde = (W_tilde + W_tilde') / 2;
    [V, D] = eig(full(W_tilde), 'vector');
    eigenvalues = max(real(D), 0);
    [eigenvalues, order] = sort(eigenvalues, 'descend');
    V = V(:, order);
    lambda1 = eigenvalues(1);
    lambda2 = eigenvalues(min(2, numel(eigenvalues)));
    rank_fraction = lambda1 / max(sum(eigenvalues), realmin);
    tail_fraction = max(0, 1 - rank_fraction);
    rank_ratio = lambda1 / max(lambda2, 1e-10);
    second_to_first_ratio = lambda2 / max(lambda1, realmin);
    relative_frobenius_rank1_error = sqrt(sum(eigenvalues(2:end).^2)) / ...
        max(norm(eigenvalues), realmin);
    numerical_rank = sum(eigenvalues > ...
        config.rank_tolerance * max(lambda1, realmin));

    [relaxed_ac, relaxed_nsre] = local_relaxed_scores( ...
        W_tilde, HB, HD, desired);

    principal = V(1:end-1, 1) / max_complex(V(end, 1));
    principal = local_limit_power(principal, config.power_bound);
    [principal_ac, principal_nsre] = local_scores( ...
        principal, HB, HD, desired);

    gaussian_ac = NaN;
    gaussian_nsre = NaN;
    best_random = nan(size(principal));
    random_ac = nan(config.randomization_trials, 1);
    random_nsre = nan(config.randomization_trials, 1);
    run_randomization = config.randomize_rank_one_solutions || ...
        rank_fraction < config.rank_fraction_threshold;
    if run_randomization
        lifted_mean = W_tilde(1:end-1, end);
        lifted_covariance = W_tilde(1:end-1, 1:end-1) - ...
            lifted_mean * lifted_mean';
        lifted_covariance = ...
            (lifted_covariance + lifted_covariance') / 2;
        [covariance_vectors, covariance_values] = eig( ...
            full(lifted_covariance), 'vector');
        covariance_values = max(real(covariance_values), 0);
        covariance_square_root = covariance_vectors * ...
            diag(sqrt(covariance_values));
        best_objective = local_nominal_proxy(principal, scaled_problem);
        best_random = principal;
        for trial = 1:config.randomization_trials
            random_vector = (randn(size(lifted_mean)) + ...
                1i * randn(size(lifted_mean))) / sqrt(2);
            candidate = lifted_mean + ...
                covariance_square_root * random_vector;
            candidate = local_limit_power(candidate, config.power_bound);
            objective = local_nominal_proxy(candidate, scaled_problem);
            [random_ac(trial), random_nsre(trial)] = local_scores( ...
                candidate, HB, HD, desired);
            if objective < best_objective
                best_objective = objective;
                best_random = candidate;
            end
        end
        [gaussian_ac, gaussian_nsre] = local_scores( ...
            best_random, HB, HD, desired);
    end
    if rank_fraction >= config.rank_fraction_threshold
        production = principal;
        production_mode = "Principal eigenvector";
    else
        production = best_random;
        production_mode = "Gaussian randomization";
    end
    [production_ac, production_nsre] = local_scores( ...
        production, HB, HD, desired);

    recovery = struct('eigenvalues', eigenvalues, ...
        'rank_fraction', rank_fraction, 'tail_fraction', tail_fraction, ...
        'rank_ratio', rank_ratio, ...
        'second_to_first_ratio', second_to_first_ratio, ...
        'relative_frobenius_rank1_error', ...
        relative_frobenius_rank1_error, ...
        'numerical_rank', numerical_rank, ...
        'production_mode', production_mode, ...
        'randomization_executed', run_randomization, ...
        'relaxed_ac_db', relaxed_ac, 'relaxed_nsre_db', relaxed_nsre, ...
        'production_filter', production, ...
        'production_ac_db', production_ac, ...
        'production_nsre_db', production_nsre, ...
        'production_minus_relaxed_ac_db', production_ac - relaxed_ac, ...
        'production_minus_relaxed_nsre_db', production_nsre - relaxed_nsre, ...
        'principal_filter', principal, 'principal_ac_db', principal_ac, ...
        'principal_nsre_db', principal_nsre, ...
        'gaussian_filter', best_random, 'gaussian_ac_db', gaussian_ac, ...
        'gaussian_nsre_db', gaussian_nsre, ...
        'gaussian_minus_principal_ac_db', gaussian_ac - principal_ac, ...
        'gaussian_minus_principal_nsre_db', gaussian_nsre - principal_nsre, ...
        'random_ac_p10_db', prctile(random_ac, 10), ...
        'random_ac_p90_db', prctile(random_ac, 90), ...
        'random_nsre_p10_db', prctile(random_nsre, 10), ...
        'random_nsre_p90_db', prctile(random_nsre, 90));
end

function [ac_db, nsre_db] = local_relaxed_scores(W_tilde, HB, HD, desired)
    L = size(HB, 2);
    W = W_tilde(1:L, 1:L);
    w = W_tilde(1:L, L+1);
    bright_energy = max(real(trace(HB * W * HB')), 0);
    dark_energy = max(real(trace(HD * W * HD')), 0);
    ratio = size(HD, 1) * bright_energy / max( ...
        size(HB, 1) * dark_energy, realmin);
    ac_db = 10 * log10(max(real(ratio), realmin));
    error_energy = real(bright_energy - ...
        2 * desired' * HB * w + desired' * desired);
    error_energy = max(error_energy, 0);
    nsre_db = 10 * log10(max(error_energy / ...
        max(norm(desired)^2, realmin), realmin));
end

function value = max_complex(value)
    if abs(value) < realmin
        value = realmin;
    end
end

function w = local_limit_power(w, power_bound)
    if norm(w)^2 > power_bound
        w = w * sqrt(power_bound) / norm(w);
    end
end

function objective = local_nominal_proxy(w, problem)
    bright = problem.HB * w;
    error = problem.desired - bright;
    dark = problem.HD * w;
    objective = real(error' * error - problem.rho * (bright' * bright) ...
        + problem.gamma * (dark' * dark));
end

function [ac_db, nsre_db] = local_scores(w, HB, HD, desired)
    ac_db = local_ac_db(w, HB, HD);
    nsre_db = 10 * log10(max(norm(HB*w-desired)^2 / ...
        max(norm(desired)^2, realmin), realmin));
end

function ac_db = local_ac_db(w, HB, HD)
    ratio = size(HD, 1) * norm(HB*w)^2 / max( ...
        size(HB, 1) * norm(HD*w)^2, realmin);
    ac_db = 10 * log10(max(real(ratio), realmin));
end

function solved = local_status_is_solved(status)
    solved = strcmp(status, 'Solved') || ...
        strcmp(status, 'Inaccurate/Solved');
end

function row = local_unsolved_row(frequency_hz, epsilon_B, epsilon_D, ...
        alpha, output)
    values = num2cell(nan(1, numel(local_variable_names())));
    values{1} = frequency_hz;
    values{3} = epsilon_B;
    values{4} = epsilon_D;
    values{5} = alpha;
    values{6} = string(output.status);
    values{7} = output.solve_time_seconds;
    values{15} = "None";
    row = cell2table(values, 'VariableNames', local_variable_names());
end

function names = local_variable_names()
    names = {'FrequencyHz', 'RelativeUncertaintyLevel', 'EpsilonB', ...
        'EpsilonD', 'Alpha', 'Status', 'SolveTimeSeconds', ...
        'NormalizationMultiplier', 'RankFraction', 'TailFraction', ...
        'RankRatio', 'SecondToFirstRatio', ...
        'RelativeFrobeniusRank1Error', 'NumericalRank', ...
        'ProductionRecoveryMode', 'RelaxedACdB', 'RelaxedNSREdB', ...
        'ProductionACdB', 'ProductionNSREdB', ...
        'ProductionMinusRelaxedACdB', ...
        'ProductionMinusRelaxedNSREdB', 'PrincipalACdB', ...
        'PrincipalNSREdB', 'GaussianACdB', 'GaussianNSREdB', ...
        'GaussianMinusPrincipalACdB', ...
        'GaussianMinusPrincipalNSREdB', 'RandomACP10dB', ...
        'RandomACP90dB', 'RandomNSREP10dB', 'RandomNSREP90dB'};
end

function band_summary = local_band_summary(summary, config)
    definitions = { ...
        "Full_100_4000", config.frequency_range_hz(1), ...
            config.frequency_range_hz(2); ...
        "Primary_100_1600", 100, 1600; ...
        "Stress_1650_4000", 1650, 4000};
    rows = cell(size(definitions, 1), 1);
    for band_no = 1:size(definitions, 1)
        name = definitions{band_no, 1};
        lower_hz = definitions{band_no, 2};
        upper_hz = definitions{band_no, 3};
        selected = summary.FrequencyHz >= lower_hz & ...
            summary.FrequencyHz <= upper_hz;
        band = summary(selected, :);
        solved = band.Status == "Solved" | ...
            band.Status == "Inaccurate/Solved";
        band = band(solved, :);
        ac_gap = abs(band.ProductionMinusRelaxedACdB);
        nsre_gap = abs(band.ProductionMinusRelaxedNSREdB);
        rows{band_no} = table(name, lower_hz, upper_hz, ...
            sum(selected), height(band), ...
            local_min(band.RankFraction), ...
            local_median(band.RankFraction), ...
            local_max(band.TailFraction), ...
            local_min(band.RankRatio), ...
            sum(band.ProductionRecoveryMode == "Gaussian randomization"), ...
            local_median(ac_gap), local_max(ac_gap), ...
            local_median(nsre_gap), local_max(nsre_gap), ...
            'VariableNames', {'Band', 'LowerHz', 'UpperHz', ...
            'RequestedFrequencyCount', 'SolvedFrequencyCount', ...
            'MinimumRankFraction', 'MedianRankFraction', ...
            'MaximumTailFraction', 'MinimumRankRatio', ...
            'GaussianRecoveryCount', 'MedianAbsoluteACGapdB', ...
            'MaximumAbsoluteACGapdB', 'MedianAbsoluteNSREGapdB', ...
            'MaximumAbsoluteNSREGapdB'});
    end
    band_summary = vertcat(rows{:});
end

function value = local_min(values)
    values = values(isfinite(values));
    if isempty(values), value = NaN; else, value = min(values); end
end

function value = local_median(values)
    values = values(isfinite(values));
    if isempty(values), value = NaN; else, value = median(values); end
end

function value = local_max(values)
    values = values(isfinite(values));
    if isempty(values), value = NaN; else, value = max(values); end
end
