function results = run_wcracc_cross_term_ablation(user_config)
%RUN_WCRACC_CROSS_TERM_ABLATION Controlled Full-WCRACC cross-term ablation.
%
% This is the reviewer-facing experiment entry point. It compares three
% robust ACC models under one common ATF-level uncertainty calibration:
%
%   1) ConventionalWCRACC:
%      correlation-matrix uncertainty with gamma_Z = eta_Z^2;
%   2) QOnlyNoCT:
%      strict ATF expansion retaining DeltaH_Z^H*DeltaH_Z only;
%   3) FullWCRACC:
%      complete ATF expansion, solved by the decomposed SDP and followed
%      by rank-one recovery.
%
% The design radii depend only on the nominal 22.5 deg C ATFs:
%
%   eta_Z(f) = nu * ||H_Z(f,22.5)||_F.
%
% The remaining 50 temperature ATFs are not used to select the uncertainty
% radius. All 51 ATFs from 20 to 25 deg C are used only for out-of-design
% evaluation. This separation avoids defining the design uncertainty from
% the same test cases used to claim robustness.
%
% The agreed primary relative uncertainty level is nu=0.005:
%   results = run_wcracc_cross_term_ablation();
%
% Useful lightweight checks:
%   cfg = struct('relative_uncertainty_levels', selected_levels, ...
%                'requested_frequencies_hz', 1000, ...
%                'run_full_sdp', false);
%   results = run_wcracc_cross_term_ablation(cfg);
%
% Relative-uncertainty sensitivity example:
%   cfg = struct('relative_uncertainty_levels', ...
%       [1e-3, 2e-3, 3e-3, 5e-3, 7.5e-3, 1e-2, 2e-2, 5e-2]);
%   results = run_wcracc_cross_term_ablation(cfg);
%
% Full frequency-grid run (computationally expensive):
%   cfg = struct('requested_frequencies_hz', []);
%   results = run_wcracc_cross_term_ablation(cfg);
%
% The normalized one-SDP formulation is attempted first. If a large,
% ill-conditioned case is reported as Unbounded or otherwise unsolved,
% config.fallback_to_bisection=true reruns that uncertainty-level batch
% with trace(W)=1
% and fixed-alpha bisection. The fallback is slower but numerically safer.
%
% Outputs are written to a timestamped directory under
% RQ_response/wcracc_cross_term_ablation_results. The principal files are:
%   experiment_results.mat       complete configuration and solver output;
%   design_diagnostics.csv       model objectives, SDR bounds, and ranks;
%   temperature_ac_metrics.csv   AC at every frequency and temperature;
%   temperature_ac_summary.csv   mean/minimum/range across temperature;
%   ac_summary_relative_level_*.png compact comparison plots.
%
% Because R_B can be severely ill-conditioned, the diagnostics report both
% the strict raw-matrix relative-level limit and a regularized effective
% limit based
% on R_B + delta_B*I, where
% delta_Z = diagonal_loading_ratio*lambda_max(R_Z). The same zone-wise
% loading is included in Conventional, Q-only, and Full-WCRACC so the
% cross-term comparison remains matched.

    if nargin < 1 || isempty(user_config)
        user_config = struct();
    end

    response_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(response_dir);
    config = local_default_config(project_root, response_dir);
    config = local_apply_overrides(config, user_config);
    local_validate_config(config);

    % Add explicit source folders instead of genpath(src). The repository
    % contains two evaluate_performance.m files with incompatible APIs.
    addpath(fullfile(project_root, 'src', 'algorithm'));
    addpath(fullfile(project_root, 'src', 'utils'));
    addpath(fullfile(project_root, 'src', 'evaluations'), '-begin');

    if config.run_full_sdp && exist('cvx_begin', 'file') == 0
        error('WCRACCAblation:CVXMissing', ...
            'CVX is required when config.run_full_sdp is true.');
    end

    rng(config.random_seed, 'twister');
    output_dir = local_prepare_output_directory(config, response_dir);
    fprintf('WCRACC cross-term ablation output: %s\n', output_dir);

    %% Step 1: Load the frequency grid and the nominal design ATFs.
    % Only the 22.5 deg C control-point ATFs are used to design filters and
    % calculate eta_Z. No off-nominal ATF enters this step.
    data_dir = fullfile(project_root, 'data', 'SimulateRIR', 'temperature');
    para_file = fullfile(data_dir, 'para.mat');
    para_data = load(para_file, 'para');
    all_frequencies = para_data.para.freq_params.target_freqs(:).';
    frequency_indices = local_frequency_indices( ...
        all_frequencies, config.requested_frequencies_hz);
    frequencies_hz = all_frequencies(frequency_indices);

    nominal_file = get_data_filename( ...
        data_dir, 'temperature', config.nominal_temperature_c);
    nominal_data = load(nominal_file, 'ATF_BZ', 'ATF_DZ');
    HB_nominal = nominal_data.ATF_BZ.ctrl(:, :, frequency_indices);
    HD_nominal = nominal_data.ATF_DZ.ctrl(:, :, frequency_indices);
    clear nominal_data

    [MB, L, num_frequencies] = size(HB_nominal);
    MD = size(HD_nominal, 1);
    fprintf('Nominal problem: MB=%d, MD=%d, L=%d, frequencies=%d.\n', ...
        MB, MD, L, num_frequencies);

    %% Step 2: Calculate admissible relative-uncertainty diagnostics.
    % relative_level_all_positive_B guarantees strictly positive worst-case
    % bright-zone energy in every loudspeaker direction. The corresponding
    % any-positive limit only guarantees that at least one such direction
    % exists. These are diagnostics and are not silently imposed by the code.
    [relative_level_all_positive_B, ...
        relative_level_all_positive_B_regularized, ...
        relative_level_any_positive_B, bright_regularization_delta] = ...
        local_relative_level_limits(HB_nominal, ...
        config.diagonal_loading_ratio);
    dark_regularization_delta = local_loading_deltas( ...
        HD_nominal, config.diagonal_loading_ratio);

    num_levels = numel(config.relative_uncertainty_levels);
    method_names = {'ConventionalWCRACC', 'QOnlyNoCT'};
    if config.run_full_sdp
        method_names{end + 1} = 'FullWCRACC';
    end
    num_methods = numel(method_names);
    num_temperatures = numel(config.test_temperatures_c);

    filter_sets = cell(num_levels, 1);
    solver_outputs = cell(num_levels, 1);
    eta_sets = cell(num_levels, 1);

    diagnostic_rows = num_levels * num_frequencies * num_methods;
    diagnostics = local_empty_diagnostics(diagnostic_rows);
    diagnostic_row = 0;

    %% Step 3: Design all filters for each predeclared uncertainty level.
    % Conventional and Q-only have generalized-eigenvalue closed forms.
    % Full-WCRACC calls the existing decomposed SDP implementation, which
    % uses the s_D=1 homogeneous normalization and therefore solves one SDP
    % per frequency rather than performing an outer alpha bisection.
    % Q-only is solved exactly under its own model, whereas Full-WCRACC is
    % compared using its recovered rank-one filter, never its relaxed SDP
    % upper bound. Consequently, a measured Full-WCRACC improvement cannot
    % be attributed to reporting the relaxation value as realizable AC.
    for level_index = 1:num_levels
        relative_level = config.relative_uncertainty_levels(level_index);
        [etaB, etaD] = local_relative_radii( ...
            HB_nominal, HD_nominal, relative_level);
        eta_sets{level_index} = struct('B', etaB, 'D', etaD);

        if any(relative_level >= relative_level_all_positive_B)
            warning('WCRACCAblation:RelativeLevelOutsideStrictRange', ...
                ['nu=%.4g exceeds the strict raw-R_B all-direction ', ...
                 'limit at one or more frequencies. The regularized ', ...
                 'effective limit is reported separately.'], relative_level);
        end
        if config.diagonal_loading_ratio > 0 && ...
                any(relative_level >= ...
                relative_level_all_positive_B_regularized)
            warning('WCRACCAblation:RelativeLevelOutsideRegularizedRange', ...
                ['nu=%.4g also exceeds the regularized all-direction ', ...
                 'limit at one or more selected frequencies.'], ...
                relative_level);
        end
        if any(relative_level >= relative_level_any_positive_B)
            warning('WCRACCAblation:RelativeLevelOutsideAnyDirectionRange', ...
                ['nu=%.4g removes the positive robust bright-zone ', ...
                 'guarantee at one or more selected frequencies.'], ...
                relative_level);
        end

        matched_para = struct();
        matched_para.full_racc.eta = struct('B', etaB, 'D', etaD);
        matched_para.full_racc.diagonal_loading_ratio = ...
            config.diagonal_loading_ratio;

        fprintf('Designing relative level %.4g: Conventional WCRACC...\n', ...
            relative_level);
        conventional = ConventionalRACCMatched( ...
            HB_nominal, HD_nominal, matched_para);

        fprintf('Designing relative level %.4g: strict Q-only/NoCT...\n', ...
            relative_level);
        qonly = local_design_qonly( ...
            HB_nominal, HD_nominal, etaB, etaD, ...
            config.diagonal_loading_ratio);

        if config.run_full_sdp
            fprintf('Designing relative level %.4g: Full-WCRACC SDP...\n', ...
                relative_level);
            full_para = matched_para;
            full_para.full_racc.formulation = 'decomposed';
            full_para.full_racc.solution_method = 'normalized';
            full_para.full_racc.solver = config.cvx_solver;
            full_para.full_racc.mosek_threads = config.mosek_threads;
            full_para.full_racc.parallel_workers = config.parallel_workers;
            full_para.full_racc.normalization_multipliers = ...
                config.normalization_multipliers;
            full_para.full_racc.checkpoint_dir = fullfile(output_dir, ...
                sprintf('full_wcracc_checkpoints_nu_%s', ...
                local_number_tag(relative_level)));
            full_para.full_racc.randomization_trials = ...
                config.randomization_trials;
            full_para.full_racc.random_seed = ...
                config.random_seed + level_index;
            full_para.full_racc.refinement_iterations = ...
                config.refinement_iterations;
            full_normalized = FullCrossTermRACC( ...
                HB_nominal, HD_nominal, full_para);
            full_bisection = [];
            full_used_fallback = false;
            normalized_solved = cellfun(@local_status_is_solved, ...
                full_normalized.status);
            if any(~normalized_solved) && config.fallback_to_bisection
                warning('WCRACCAblation:NormalizedSDPFailed', ...
                    ['The normalized SDP was not solved at every selected ', ...
                     'frequency. Retrying this level with trace(W)=1 and ', ...
                     'fixed-alpha bisection.']);
                full_para.full_racc.solution_method = 'bisection';
                full_para.full_racc.bisection_tolerance_db = ...
                    config.bisection_tolerance_db;
                full_para.full_racc.max_bisection_iterations = ...
                    config.max_bisection_iterations;
                full_bisection = FullCrossTermRACC( ...
                    HB_nominal, HD_nominal, full_para);
                full = full_bisection;
                full_used_fallback = true;
            else
                full = full_normalized;
            end
        else
            full = [];
            full_normalized = [];
            full_bisection = [];
            full_used_fallback = false;
        end

        filters = struct();
        filters.ConventionalWCRACC = conventional.w;
        filters.QOnlyNoCT = qonly.w;
        if config.run_full_sdp
            filters.FullWCRACC = full.w;
        end
        filter_sets{level_index} = filters;
        solver_outputs{level_index} = struct( ...
            'conventional', conventional, 'qonly', qonly, 'full', full, ...
            'full_normalized_attempt', full_normalized, ...
            'full_bisection_fallback', full_bisection, ...
            'full_used_fallback', full_used_fallback);

        % Save one checkpoint per level so a later failure does not discard
        % already completed, computationally expensive SDP solutions.
        checkpoint_file = fullfile(output_dir, sprintf( ...
            'design_checkpoint_relative_level_%s.mat', ...
            local_number_tag(relative_level)));
        save(checkpoint_file, 'relative_level', 'etaB', 'etaD', 'filters', ...
            'conventional', 'qonly', 'full', 'full_normalized', ...
            'full_bisection', 'full_used_fallback', '-v7.3');

        for frequency_index = 1:num_frequencies
            for method_index = 1:num_methods
                diagnostic_row = diagnostic_row + 1;
                method = method_names{method_index};
                w = filters.(method)(:, frequency_index);
                exact_full_db = local_exact_full_ball_score( ...
                    w, HB_nominal(:, :, frequency_index), ...
                    HD_nominal(:, :, frequency_index), ...
                    etaB(frequency_index), etaD(frequency_index), ...
                    bright_regularization_delta(frequency_index), ...
                    dark_regularization_delta(frequency_index));

                sdp_upper_db = NaN;
                rank_fraction = NaN;
                sdr_gap_db = NaN;
                status = "Closed form";
                if strcmp(method, 'ConventionalWCRACC')
                    internal_db = conventional.correlation_model_ac_db( ...
                        frequency_index);
                elseif strcmp(method, 'QOnlyNoCT')
                    internal_db = qonly.model_ac_db(frequency_index);
                else
                    internal_db = full.recovered_ac_db(frequency_index);
                    sdp_upper_db = full.relaxation_ac_db(frequency_index);
                    rank_fraction = full.rank_fraction(frequency_index);
                    sdr_gap_db = sdp_upper_db - exact_full_db;
                    status = string(full.status{frequency_index});
                    if full_used_fallback
                        status = "Bisection fallback: " + status;
                    end
                end

                diagnostics.RelativeUncertaintyLevel(diagnostic_row) = ...
                    relative_level;
                diagnostics.FrequencyHz(diagnostic_row) = ...
                    frequencies_hz(frequency_index);
                diagnostics.Algorithm(diagnostic_row) = string(method);
                diagnostics.EtaB(diagnostic_row) = etaB(frequency_index);
                diagnostics.EtaD(diagnostic_row) = etaD(frequency_index);
                diagnostics.GammaB(diagnostic_row) = etaB(frequency_index)^2;
                diagnostics.GammaD(diagnostic_row) = etaD(frequency_index)^2;
                diagnostics.RelativeLevelAllPositiveB(diagnostic_row) = ...
                    relative_level_all_positive_B(frequency_index);
                diagnostics.BrightRegularizationDelta(diagnostic_row) = ...
                    bright_regularization_delta(frequency_index);
                diagnostics.DarkRegularizationDelta(diagnostic_row) = ...
                    dark_regularization_delta(frequency_index);
                diagnostics.RelativeLevelAllPositiveBRegularized( ...
                    diagnostic_row) = ...
                    relative_level_all_positive_B_regularized( ...
                    frequency_index);
                diagnostics.RelativeLevelAnyPositiveB(diagnostic_row) = ...
                    relative_level_any_positive_B(frequency_index);
                diagnostics.InternalObjectiveACdB(diagnostic_row) = internal_db;
                diagnostics.ExactRegularizedATFBallACdB( ...
                    diagnostic_row) = exact_full_db;
                diagnostics.SDPUpperBoundACdB(diagnostic_row) = sdp_upper_db;
                diagnostics.RankFraction(diagnostic_row) = rank_fraction;
                diagnostics.SDRGapdB(diagnostic_row) = sdr_gap_db;
                diagnostics.Status(diagnostic_row) = status;
            end
        end
    end

    design_diagnostics = struct2table(diagnostics);
    writetable(design_diagnostics, fullfile( ...
        output_dir, 'design_diagnostics.csv'));

    %% Step 4: Evaluate unchanged filters on all 51 temperature ATFs.
    % Both control-point AC and independent evaluation-point AC are stored.
    % RelativeDeltaH is reported only to describe the physical mismatch; it
    % is never fed back into eta_Z or used to redesign a filter.
    metric_rows = num_levels * num_temperatures * ...
        num_frequencies * num_methods;
    metrics = local_empty_metrics(metric_rows);
    metric_row = 0;

    for temperature_index = 1:num_temperatures
        temperature_c = config.test_temperatures_c(temperature_index);
        test_file = get_data_filename(data_dir, 'temperature', temperature_c);
        test_data = load(test_file, 'ATF_BZ', 'ATF_DZ');
        HB_control = test_data.ATF_BZ.ctrl(:, :, frequency_indices);
        HD_control = test_data.ATF_DZ.ctrl(:, :, frequency_indices);
        HB_evaluation = test_data.ATF_BZ.eval(:, :, frequency_indices);
        HD_evaluation = test_data.ATF_DZ.eval(:, :, frequency_indices);
        clear test_data

        [relative_delta_B, relative_delta_D] = local_relative_perturbations( ...
            HB_control, HD_control, HB_nominal, HD_nominal);

        for level_index = 1:num_levels
            filters = filter_sets{level_index};
            for method_index = 1:num_methods
                method = method_names{method_index};
                w_all = filters.(method);
                for frequency_index = 1:num_frequencies
                    metric_row = metric_row + 1;
                    w = w_all(:, frequency_index);
                    metrics.RelativeUncertaintyLevel(metric_row) = ...
                        config.relative_uncertainty_levels(level_index);
                    metrics.TemperatureC(metric_row) = temperature_c;
                    metrics.FrequencyHz(metric_row) = ...
                        frequencies_hz(frequency_index);
                    metrics.Algorithm(metric_row) = string(method);
                    metrics.ControlACdB(metric_row) = local_ac_db( ...
                        w, HB_control(:, :, frequency_index), ...
                        HD_control(:, :, frequency_index));
                    metrics.EvaluationACdB(metric_row) = local_ac_db( ...
                        w, HB_evaluation(:, :, frequency_index), ...
                        HD_evaluation(:, :, frequency_index));
                    metrics.RelativeDeltaHB(metric_row) = ...
                        relative_delta_B(frequency_index);
                    metrics.RelativeDeltaHD(metric_row) = ...
                        relative_delta_D(frequency_index);
                end
            end
        end
        fprintf('Evaluated temperature %.1f deg C (%d/%d).\n', ...
            temperature_c, temperature_index, num_temperatures);
    end

    temperature_metrics = struct2table(metrics);
    writetable(temperature_metrics, fullfile( ...
        output_dir, 'temperature_ac_metrics.csv'));

    %% Step 5: Summarize, plot, and save all reproducibility information.
    temperature_summary = local_temperature_summary( ...
        temperature_metrics, config.relative_uncertainty_levels, ...
        frequencies_hz, method_names);
    writetable(temperature_summary, fullfile( ...
        output_dir, 'temperature_ac_summary.csv'));
    local_plot_summaries(temperature_summary, ...
        config.relative_uncertainty_levels, ...
        method_names, output_dir);

    results = struct();
    results.config = config;
    results.output_dir = output_dir;
    results.frequency_indices = frequency_indices;
    results.frequencies_hz = frequencies_hz;
    results.relative_level_limits = struct( ...
        'all_positive_B', relative_level_all_positive_B, ...
        'all_positive_B_regularized', ...
        relative_level_all_positive_B_regularized, ...
        'any_positive_B', relative_level_any_positive_B, ...
        'regularization_delta_B', bright_regularization_delta, ...
        'regularization_delta_D', dark_regularization_delta, ...
        'regularization_ratio', config.diagonal_loading_ratio);
    results.eta = eta_sets;
    results.filters = filter_sets;
    results.solver_outputs = solver_outputs;
    results.design_diagnostics = design_diagnostics;
    results.temperature_metrics = temperature_metrics;
    results.temperature_summary = temperature_summary;
    save(fullfile(output_dir, 'experiment_results.mat'), 'results', '-v7.3');

    fprintf('Completed WCRACC cross-term ablation: %s\n', output_dir);
end

function config = local_default_config(project_root, response_dir)
% Keep reviewer-facing choices together. requested_frequencies_hz=[] selects
% the complete stored grid; the four-frequency default is safer for review.
    config = struct();
    config.project_root = project_root;
    config.nominal_temperature_c = 22.5;
    config.test_temperatures_c = 20.0:0.1:25.0;
    config.requested_frequencies_hz = [500, 1000, 2000, 4000];
    % Primary reviewer-facing setting: a 0.5% relative ATF radius.
    config.relative_uncertainty_levels = 5e-3;
    config.diagonal_loading_ratio = 1e-6;
    config.run_full_sdp = true;
    config.fallback_to_bisection = true;
    config.bisection_tolerance_db = 0.5;
    config.max_bisection_iterations = 10;
    config.cvx_solver = 'mosek';
    config.mosek_threads = 4;
    config.parallel_workers = 1;
    config.normalization_multipliers = [1, 5, 0.1, 10, 0.01];
    config.randomization_trials = 300;
    config.refinement_iterations = 300;
    config.random_seed = 20260822;
    config.output_dir = '';
    config.output_root = fullfile( ...
        response_dir, 'wcracc_cross_term_ablation_results');
end

function config = local_apply_overrides(config, overrides)
    if ~isstruct(overrides) || ~isscalar(overrides)
        error('WCRACCAblation:InvalidConfig', ...
            'user_config must be a scalar structure.');
    end
    names = fieldnames(overrides);
    for k = 1:numel(names)
        if ~isfield(config, names{k})
            error('WCRACCAblation:UnknownConfigField', ...
                'Unknown configuration field: %s', names{k});
        end
        config.(names{k}) = overrides.(names{k});
    end
end

function local_validate_config(config)
    validateattributes(config.nominal_temperature_c, {'numeric'}, ...
        {'scalar', 'real', 'finite'});
    validateattributes(config.test_temperatures_c, {'numeric'}, ...
        {'vector', 'real', 'finite', 'nonempty'});
    validateattributes(config.relative_uncertainty_levels, {'numeric'}, ...
        {'vector', 'real', 'finite', 'positive', 'nonempty'});
    validateattributes(config.diagonal_loading_ratio, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(config.mosek_threads, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(config.parallel_workers, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(config.normalization_multipliers, {'numeric'}, ...
        {'vector', 'real', 'finite', 'positive', 'nonempty'});
    validateattributes(config.randomization_trials, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative'});
    validateattributes(config.run_full_sdp, {'logical', 'numeric'}, ...
        {'scalar'});
    validateattributes(config.fallback_to_bisection, {'logical', 'numeric'}, ...
        {'scalar'});
    validateattributes(config.bisection_tolerance_db, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(config.max_bisection_iterations, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    if ~isempty(config.requested_frequencies_hz)
        validateattributes(config.requested_frequencies_hz, {'numeric'}, ...
            {'vector', 'real', 'finite', 'positive'});
    end
end

function output_dir = local_prepare_output_directory(config, response_dir)
    if isempty(config.output_dir)
        timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        output_dir = fullfile(config.output_root, ['run_' timestamp]);
    else
        output_dir = char(string(config.output_dir));
        if ~local_is_absolute_path(output_dir)
            output_dir = fullfile(response_dir, output_dir);
        end
    end
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
end

function solved = local_status_is_solved(status)
    status = char(string(status));
    solved = strcmp(status, 'Solved') || strcmp(status, 'Inaccurate/Solved');
end

function tf = local_is_absolute_path(path_value)
    tf = ~isempty(regexp(path_value, '^[A-Za-z]:[\\/]', 'once')) || ...
        startsWith(path_value, '\\') || startsWith(path_value, '/');
end

function indices = local_frequency_indices(all_frequencies, requested)
    if isempty(requested)
        indices = 1:numel(all_frequencies);
        return;
    end
    indices = zeros(size(requested));
    for k = 1:numel(requested)
        [distance, indices(k)] = min(abs(all_frequencies - requested(k)));
        if distance > 1e-8
            warning('WCRACCAblation:NearestFrequencyUsed', ...
                'Requested %.3f Hz; using stored frequency %.3f Hz.', ...
                requested(k), all_frequencies(indices(k)));
        end
    end
    indices = unique(indices, 'stable');
end

function [strict_limit, regularized_limit, upper_limit, deltaB] = ...
        local_relative_level_limits(HB, regularization_ratio)
    num_frequencies = size(HB, 3);
    strict_limit = nan(1, num_frequencies);
    regularized_limit = nan(1, num_frequencies);
    upper_limit = nan(1, num_frequencies);
    deltaB = nan(1, num_frequencies);
    for k = 1:num_frequencies
        H = HB(:, :, k);
        singular_values = svd(H, 'econ');
        trace_RB = norm(H, 'fro')^2;
        lambda_max_RB = max(singular_values)^2;
        deltaB(k) = regularization_ratio * lambda_max_RB;
        strict_limit(k) = min(singular_values) / ...
            max(sqrt(trace_RB), realmin);
        regularized_limit(k) = sqrt( ...
            (min(singular_values)^2 + deltaB(k)) / ...
            max(trace_RB, realmin));
        upper_limit(k) = max(singular_values) / ...
            max(sqrt(trace_RB), realmin);
    end
end

function delta = local_loading_deltas(H, regularization_ratio)
    num_frequencies = size(H, 3);
    delta = nan(1, num_frequencies);
    for k = 1:num_frequencies
        delta(k) = regularization_ratio * norm(H(:, :, k), 2)^2;
    end
end

function [etaB, etaD] = local_relative_radii(HB, HD, relative_level)
    num_frequencies = size(HB, 3);
    etaB = zeros(1, num_frequencies);
    etaD = zeros(1, num_frequencies);
    for k = 1:num_frequencies
        etaB(k) = relative_level * norm(HB(:, :, k), 'fro');
        etaD(k) = relative_level * norm(HD(:, :, k), 'fro');
    end
end

function result = local_design_qonly( ...
        HB, HD, etaB, etaD, diagonal_loading_ratio)
% Strict NoCT model:
%   Rtilde_Z = R_Z + DeltaH_Z^H*DeltaH_Z.
% The bright-zone minimum is attained at DeltaH_B=0, so etaB does not
% appear in the generalized eigenproblem. The dark-zone maximum contributes
% etaD^2*I. etaB remains in the common full-ball post-evaluation.
    [MB, L, num_frequencies] = size(HB);
    MD = size(HD, 1);
    w_all = zeros(L, num_frequencies);
    alpha_all = nan(1, num_frequencies);
    ac_db_all = nan(1, num_frequencies);
    deltaB_all = nan(1, num_frequencies);
    deltaD_all = nan(1, num_frequencies);

    for k = 1:num_frequencies
        HBi = HB(:, :, k);
        HDi = HD(:, :, k);
        scale = max(norm([HBi; HDi], 'fro'), realmin);
        HBs = HBi / scale;
        HDs = HDi / scale;
        etaDs = etaD(k) / scale;
        deltaB = diagonal_loading_ratio * norm(HBi, 2)^2;
        deltaD = diagonal_loading_ratio * norm(HDi, 2)^2;
        deltaBs = deltaB / scale^2;
        deltaDs = deltaD / scale^2;

        RB = HBs' * HBs;
        RD = HDs' * HDs;
        RB = (RB + RB') / 2;
        RD = (RD + RD') / 2;
        A = RB + deltaBs * eye(L);
        B = RD + (deltaDs + etaDs^2) * eye(L);
        numerical_floor = max(1e-12 * real(trace(B)) / L, eps);
        [V, D] = eig(A, B + numerical_floor * eye(L));
        [~, index] = max(real(diag(D)));
        w = V(:, index);
        w = w / max(norm(w), realmin);

        numerator = real(w' * A * w);
        denominator = real(w' * B * w);
        alpha = (MD / MB) * numerator / max(denominator, realmin);
        w_all(:, k) = w;
        alpha_all(k) = alpha;
        ac_db_all(k) = 10 * log10(max(alpha, realmin));
        deltaB_all(k) = deltaB;
        deltaD_all(k) = deltaD;
    end

    result = struct('w', w_all, 'etaB', etaB, 'etaD', etaD, ...
        'model_alpha', alpha_all, 'model_ac_db', ac_db_all, ...
        'diagonal_loading', struct('B', deltaB_all, 'D', deltaD_all, ...
        'ratio', diagonal_loading_ratio));
end

function score_db = local_exact_full_ball_score( ...
        w, HB, HD, etaB, etaD, deltaB, deltaD)
    if any(~isfinite(w)) || norm(w) <= realmin
        score_db = NaN;
        return;
    end
    score_db = full_racc_worst_case_contrast( ...
        w, HB, HD, etaB, etaD, deltaB, deltaD);
end

function [relative_B, relative_D] = local_relative_perturbations( ...
        HB, HD, HB_nominal, HD_nominal)
    num_frequencies = size(HB, 3);
    relative_B = zeros(1, num_frequencies);
    relative_D = zeros(1, num_frequencies);
    for k = 1:num_frequencies
        relative_B(k) = norm(HB(:, :, k) - HB_nominal(:, :, k), 'fro') / ...
            max(norm(HB_nominal(:, :, k), 'fro'), realmin);
        relative_D(k) = norm(HD(:, :, k) - HD_nominal(:, :, k), 'fro') / ...
            max(norm(HD_nominal(:, :, k), 'fro'), realmin);
    end
end

function ac_db = local_ac_db(w, HB, HD)
    if any(~isfinite(w)) || norm(w) <= realmin
        ac_db = NaN;
        return;
    end
    bright_energy = norm(HB * w)^2;
    dark_energy = norm(HD * w)^2;
    MB = size(HB, 1);
    MD = size(HD, 1);
    alpha = (MD * bright_energy) / max(MB * dark_energy, realmin);
    ac_db = 10 * log10(max(real(alpha), realmin));
end

function data = local_empty_diagnostics(num_rows)
    data = struct();
    data.RelativeUncertaintyLevel = nan(num_rows, 1);
    data.FrequencyHz = nan(num_rows, 1);
    data.Algorithm = strings(num_rows, 1);
    data.EtaB = nan(num_rows, 1);
    data.EtaD = nan(num_rows, 1);
    data.GammaB = nan(num_rows, 1);
    data.GammaD = nan(num_rows, 1);
    data.RelativeLevelAllPositiveB = nan(num_rows, 1);
    data.BrightRegularizationDelta = nan(num_rows, 1);
    data.DarkRegularizationDelta = nan(num_rows, 1);
    data.RelativeLevelAllPositiveBRegularized = nan(num_rows, 1);
    data.RelativeLevelAnyPositiveB = nan(num_rows, 1);
    data.InternalObjectiveACdB = nan(num_rows, 1);
    data.ExactRegularizedATFBallACdB = nan(num_rows, 1);
    data.SDPUpperBoundACdB = nan(num_rows, 1);
    data.RankFraction = nan(num_rows, 1);
    data.SDRGapdB = nan(num_rows, 1);
    data.Status = strings(num_rows, 1);
end

function data = local_empty_metrics(num_rows)
    data = struct();
    data.RelativeUncertaintyLevel = nan(num_rows, 1);
    data.TemperatureC = nan(num_rows, 1);
    data.FrequencyHz = nan(num_rows, 1);
    data.Algorithm = strings(num_rows, 1);
    data.ControlACdB = nan(num_rows, 1);
    data.EvaluationACdB = nan(num_rows, 1);
    data.RelativeDeltaHB = nan(num_rows, 1);
    data.RelativeDeltaHD = nan(num_rows, 1);
end

function summary = local_temperature_summary( ...
        metrics, relative_levels, frequencies, method_names)
    num_rows = numel(relative_levels) * numel(frequencies) * ...
        numel(method_names);
    RelativeUncertaintyLevel = nan(num_rows, 1);
    FrequencyHz = nan(num_rows, 1);
    Algorithm = strings(num_rows, 1);
    MeanControlACdB = nan(num_rows, 1);
    MinimumControlACdB = nan(num_rows, 1);
    ControlACRangedB = nan(num_rows, 1);
    MeanEvaluationACdB = nan(num_rows, 1);
    MinimumEvaluationACdB = nan(num_rows, 1);
    EvaluationACRangedB = nan(num_rows, 1);
    row = 0;

    for relative_level = relative_levels(:).'
        for frequency = frequencies(:).'
            for method_index = 1:numel(method_names)
                row = row + 1;
                method = string(method_names{method_index});
                selected = abs(metrics.RelativeUncertaintyLevel - ...
                    relative_level) <= ...
                    10 * eps(max(abs(relative_level), 1)) & ...
                    metrics.FrequencyHz == frequency & ...
                    metrics.Algorithm == method;
                control = metrics.ControlACdB(selected);
                evaluation = metrics.EvaluationACdB(selected);
                RelativeUncertaintyLevel(row) = relative_level;
                FrequencyHz(row) = frequency;
                Algorithm(row) = method;
                MeanControlACdB(row) = mean(control, 'omitnan');
                MinimumControlACdB(row) = min(control, [], 'omitnan');
                ControlACRangedB(row) = ...
                    max(control, [], 'omitnan') - min(control, [], 'omitnan');
                MeanEvaluationACdB(row) = mean(evaluation, 'omitnan');
                MinimumEvaluationACdB(row) = min(evaluation, [], 'omitnan');
                EvaluationACRangedB(row) = ...
                    max(evaluation, [], 'omitnan') - ...
                    min(evaluation, [], 'omitnan');
            end
        end
    end

    summary = table(RelativeUncertaintyLevel, FrequencyHz, Algorithm, ...
        MeanControlACdB, ...
        MinimumControlACdB, ControlACRangedB, MeanEvaluationACdB, ...
        MinimumEvaluationACdB, EvaluationACRangedB);
end

function local_plot_summaries(summary, relative_levels, method_names, output_dir)
    colors = lines(numel(method_names));
    for level_index = 1:numel(relative_levels)
        relative_level = relative_levels(level_index);
        selected_level = abs(summary.RelativeUncertaintyLevel - ...
            relative_level) <= 10 * eps(max(abs(relative_level), 1));
        figure_handle = figure('Visible', 'off', 'Color', 'w', ...
            'Position', [100, 100, 1050, 430]);
        layout = tiledlayout(1, 2, 'TileSpacing', 'compact', ...
            'Padding', 'compact');

        nexttile(layout, 1);
        hold on;
        for method_index = 1:numel(method_names)
            method = string(method_names{method_index});
            selected = selected_level & summary.Algorithm == method;
            plot(summary.FrequencyHz(selected), ...
                summary.MeanEvaluationACdB(selected), '-o', ...
                'Color', colors(method_index, :), 'LineWidth', 1.4, ...
                'DisplayName', char(method));
        end
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Mean evaluation AC (dB)');
        title(sprintf('Mean across 20-25 deg C, nu=%.4g', ...
            relative_level));

        nexttile(layout, 2);
        hold on;
        for method_index = 1:numel(method_names)
            method = string(method_names{method_index});
            selected = selected_level & summary.Algorithm == method;
            plot(summary.FrequencyHz(selected), ...
                summary.MinimumEvaluationACdB(selected), '-o', ...
                'Color', colors(method_index, :), 'LineWidth', 1.4, ...
                'DisplayName', char(method));
        end
        grid on;
        xlabel('Frequency (Hz)');
        ylabel('Minimum evaluation AC (dB)');
        title('Worst observed temperature');
        legend('Location', 'best');

        exportgraphics(figure_handle, fullfile(output_dir, sprintf( ...
            'ac_summary_relative_level_%s.png', ...
            local_number_tag(relative_level))), ...
            'Resolution', 200);
        close(figure_handle);
    end
end

function tag = local_number_tag(value)
    tag = strrep(sprintf('%.6g', value), '.', 'p');
    tag = strrep(tag, '-', 'm');
    tag = strrep(tag, '+', 'p');
end
