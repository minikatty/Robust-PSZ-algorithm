function results = run_racc_pm_uncertainty_sweep(user_config)
%RUN_RACC_PM_UNCERTAINTY_SWEEP Selected-frequency RACC-PM radius study.
%
% The proposed RACC-PM method is the main subject of this experiment.
% RACC_PM_Sub (the maintained no-suffix implementation) is used throughout;
% the computationally expensive GLS implementation is deliberately excluded.
%
% The nominal 22.5 deg C control ATFs are used to design every filter. Alpha
% is calculated once from a NoCT-WCRACC control-point reference and is then
% held fixed while the common relative ATF radius nu_Z is swept. Evaluation
% ATFs are not accessed until all designs have been completed.
%
% Primary methods:
%   RACC-PM  proposed hybrid method, solved by RACC_PM_Sub with MOSEK;
%   RPM      robust pressure-matching reference at the same ATF radii;
%   wcRACC   matched conventional correlation-level robust ACC;
%   ACC-Reg  fixed small-loading energy-contrast baseline.
%
% Default constant-radius candidates:
%   nu_B = nu_D = [0.001 0.0025 0.005 0.0075 0.01 0.02 0.05].
% The manuscript rule epsilon_Z=0.01*sqrt(||R_Z||_F) is additionally run as
% a predeclared frequency-dependent candidate.
%
% Example smoke test:
%   cfg = struct('requested_frequencies_hz', 1000, ...
%       'relative_uncertainty_levels', 0.005, ...
%       'include_manuscript_rule', false, ...
%       'test_temperatures_c', 22.5);
%   run_racc_pm_uncertainty_sweep(cfg);

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

    % The project contains two functions with this name. This experiment
    % uses local metric helpers, but records the active project resolution so
    % a later full-pipeline run cannot silently switch implementations.
    active_evaluator = which('evaluate_performance');
    expected_evaluator = fullfile(project_root, 'src', 'evaluations', ...
        'evaluate_performance.m');
    if ~strcmpi(active_evaluator, expected_evaluator)
        error('RACCPMRadiusSweep:EvaluatorResolutionMismatch', ...
            'Expected %s but MATLAB resolves evaluate_performance to %s.', ...
            expected_evaluator, active_evaluator);
    end
    active_racc_pm_solver = which('RACC_PM_Sub');
    expected_solver = fullfile(project_root, 'src', 'algorithm', ...
        'RACC_PM_Sub.m');
    if ~strcmpi(active_racc_pm_solver, expected_solver)
        error('RACCPMRadiusSweep:SolverResolutionMismatch', ...
            'Expected %s but MATLAB resolves RACC_PM_Sub to %s.', ...
            expected_solver, active_racc_pm_solver);
    end
    if exist('cvx_begin', 'file') == 0
        error('RACCPMRadiusSweep:CVXMissing', ...
            'CVX with MOSEK is required for RACC_PM_Sub.');
    end

    output_dir = local_prepare_output_directory(config);
    fprintf('RACC-PM uncertainty sweep output: %s\n', output_dir);
    fprintf('Active evaluator: %s\n', active_evaluator);
    fprintf('Active RACC-PM solver: %s\n', active_racc_pm_solver);

    %% Load the selected nominal problem.
    data_dir = fullfile(project_root, 'data', 'SimulateRIR', 'temperature');
    para_data = load(fullfile(data_dir, 'para.mat'), 'para');
    all_frequencies_hz = para_data.para.freq_params.target_freqs(:).';
    frequency_indices = local_frequency_indices(all_frequencies_hz, ...
        config.requested_frequencies_hz);
    frequencies_hz = all_frequencies_hz(frequency_indices);

    nominal_file = get_data_filename(data_dir, 'temperature', ...
        config.nominal_temperature_c);
    nominal_data = load(nominal_file, 'ATF_BZ', 'ATF_DZ');
    HB_nominal = nominal_data.ATF_BZ.ctrl(:, :, frequency_indices);
    HD_nominal = nominal_data.ATF_DZ.ctrl(:, :, frequency_indices);
    clear nominal_data

    desired_file = fullfile(project_root, 'data', ...
        'ATF_desired_plane.mat');
    desired_data = load(desired_file, 'ATF_desired_plane');
    desired_pressure = desired_data.ATF_desired_plane(:, frequency_indices);
    clear desired_data

    [MB, L, num_frequencies] = size(HB_nominal);
    MD = size(HD_nominal, 1);
    fprintf('Nominal problem: MB=%d, MD=%d, L=%d, frequencies=%s Hz.\n', ...
        MB, MD, L, mat2str(frequencies_hz));

    %% Define the fixed control-derived alpha before any evaluation ATF load.
    [eta_reference_B, eta_reference_D] = local_relative_radii( ...
        HB_nominal, HD_nominal, config.alpha_reference_relative_level);
    alpha_reference_para = local_matched_para( ...
        eta_reference_B, eta_reference_D, 0);
    alpha_reference_result = NoCT_WCRACC( ...
        HB_nominal, HD_nominal, alpha_reference_para);
    alpha_reference_ac_db = nan(1, num_frequencies);
    for frequency_index = 1:num_frequencies
        alpha_reference_ac_db(frequency_index) = local_ac_db( ...
            alpha_reference_result.w(:, frequency_index), ...
            HB_nominal(:, :, frequency_index), ...
            HD_nominal(:, :, frequency_index));
    end
    alpha = 10.^((alpha_reference_ac_db - ...
        config.alpha_margin_db) / 10);
    fprintf(['Fixed alpha from nominal NoCT control points ', ...
        '(nu_ref=%.4g, margin=%.3g dB):\n'], ...
        config.alpha_reference_relative_level, config.alpha_margin_db);
    for frequency_index = 1:num_frequencies
        fprintf('  %.0f Hz: AC_ref=%.6f dB, alpha=%.9g\n', ...
            frequencies_hz(frequency_index), ...
            alpha_reference_ac_db(frequency_index), ...
            alpha(frequency_index));
    end

    %% Create the predeclared uncertainty settings.
    settings = local_build_settings(config, HB_nominal, HD_nominal);
    num_settings = numel(settings);
    fprintf('Uncertainty settings (%d): %s\n', num_settings, ...
        strjoin({settings.label}, ', '));

    %% Fixed ACC-Reg baseline.
    [acc_reg_filter, acc_reg_diagnostics] = local_design_acc_reg( ...
        HB_nominal, HD_nominal, config.acc_reg_beta, frequencies_hz);
    writetable(acc_reg_diagnostics, fullfile(output_dir, ...
        'acc_reg_design_diagnostics.csv'));

    %% Design each uncertainty-dependent method and checkpoint immediately.
    designs = cell(num_settings, 1);
    design_rows = num_settings * num_frequencies * 4;
    design_diagnostics = local_empty_design_table(design_rows);
    design_row = 0;

    for setting_index = 1:num_settings
        setting = settings(setting_index);
        checkpoint_file = fullfile(output_dir, sprintf( ...
            'design_checkpoint_%02d_%s.mat', setting_index, setting.id));

        if config.resume && isfile(checkpoint_file)
            checkpoint = load(checkpoint_file, 'design');
            design = checkpoint.design;
            fprintf('Reusing checkpoint %d/%d: %s\n', setting_index, ...
                num_settings, checkpoint_file);
        else
            fprintf('\nSetting %d/%d: %s\n', setting_index, ...
                num_settings, setting.label);

            matched_para = local_matched_para( ...
                setting.epsilon_B, setting.epsilon_D, 0);
            wcracc_result = ConventionalRACCMatched( ...
                HB_nominal, HD_nominal, matched_para);

            robust_para = struct();
            robust_para.scale = 1;
            robust_para.epsilon = struct('B', setting.epsilon_B, ...
                'D', setting.epsilon_D);
            robust_para.rho = config.rho;
            robust_para.mu = config.mu;
            robust_para.alpha = alpha;

            rpm_result = [];
            rpm_error = "";
            if config.run_rpm
                try
                    rpm_result = RPM(HB_nominal, HD_nominal, ...
                        desired_pressure, robust_para);
                catch exception
                    rpm_error = string(getReport(exception, 'extended', ...
                        'hyperlinks', 'off'));
                    warning('RACCPMRadiusSweep:RPMFailed', ...
                        'RPM failed for %s: %s', setting.label, ...
                        exception.message);
                end
            end

            racc_pm_result = [];
            racc_pm_error = "";
            try
                racc_pm_result = RACC_PM_Sub(HB_nominal, HD_nominal, ...
                    desired_pressure, robust_para, 1:num_frequencies);
            catch exception
                racc_pm_error = string(getReport(exception, 'extended', ...
                    'hyperlinks', 'off'));
                warning('RACCPMRadiusSweep:RACCPMFailed', ...
                    'RACC-PM failed for %s: %s', setting.label, ...
                    exception.message);
            end

            filters = struct();
            filters.ACC_Reg = acc_reg_filter;
            filters.wcRACC = wcracc_result.w;
            if isempty(rpm_result)
                filters.RPM = nan(L, num_frequencies);
            else
                filters.RPM = rpm_result.w;
            end
            if isempty(racc_pm_result)
                filters.RACC_PM = nan(L, num_frequencies);
            else
                filters.RACC_PM = racc_pm_result.w;
            end

            design = struct();
            design.setting = setting;
            design.filters = filters;
            design.wcracc_result = wcracc_result;
            design.rpm_result = rpm_result;
            design.rpm_error = rpm_error;
            design.racc_pm_result = racc_pm_result;
            design.racc_pm_error = racc_pm_error;
            design.alpha = alpha;
            design.alpha_reference_ac_db = alpha_reference_ac_db;
            save(checkpoint_file, 'design', '-v7.3');
        end
        designs{setting_index} = design;

        algorithm_ids = {'ACC_Reg', 'wcRACC', 'RPM', 'RACC_PM'};
        display_names = {'ACC-Reg', 'wcRACC', 'RPM', 'RACC-PM'};
        for algorithm_index = 1:numel(algorithm_ids)
            algorithm_id = algorithm_ids{algorithm_index};
            for frequency_index = 1:num_frequencies
                design_row = design_row + 1;
                w = design.filters.(algorithm_id)(:, frequency_index);
                valid_filter = all(isfinite(w)) && norm(w) > 0;
                design_diagnostics.Setting(design_row) = ...
                    string(setting.label);
                design_diagnostics.Algorithm(design_row) = ...
                    string(display_names{algorithm_index});
                design_diagnostics.FrequencyHz(design_row) = ...
                    frequencies_hz(frequency_index);
                design_diagnostics.NuB(design_row) = ...
                    setting.nu_B(frequency_index);
                design_diagnostics.NuD(design_row) = ...
                    setting.nu_D(frequency_index);
                design_diagnostics.EpsilonB(design_row) = ...
                    setting.epsilon_B(frequency_index);
                design_diagnostics.EpsilonD(design_row) = ...
                    setting.epsilon_D(frequency_index);
                design_diagnostics.Alpha(design_row) = ...
                    alpha(frequency_index);
                design_diagnostics.ValidFilter(design_row) = valid_filter;
                design_diagnostics.FilterNorm(design_row) = norm(w);
                design_diagnostics.Status(design_row) = ...
                    local_design_status(design, algorithm_id, ...
                    frequency_index, valid_filter);
                if valid_filter
                    design_diagnostics.ControlACdB(design_row) = ...
                        local_ac_db(w, HB_nominal(:, :, frequency_index), ...
                        HD_nominal(:, :, frequency_index));
                    if any(strcmp(algorithm_id, {'RPM', 'RACC_PM'}))
                        design_diagnostics.ControlNSREdB(design_row) = ...
                            local_nsre_db(w, ...
                            HB_nominal(:, :, frequency_index), ...
                            desired_pressure(:, frequency_index));
                    end
                    design_diagnostics.ControlAEdB(design_row) = ...
                        local_ae_db(w, ...
                        HB_nominal(:, :, frequency_index), ...
                        config.virtual_source_idx, config.target_spl_db);
                end
            end
        end
        writetable(design_diagnostics(1:design_row, :), fullfile( ...
            output_dir, 'design_diagnostics.csv'));
    end

    %% Held-out spatial/temperature evaluation starts here.
    algorithm_ids = {'ACC_Reg', 'wcRACC', 'RPM', 'RACC_PM'};
    display_names = {'ACC-Reg', 'wcRACC', 'RPM', 'RACC-PM'};
    num_algorithms = numel(algorithm_ids);
    num_temperatures = numel(config.test_temperatures_c);
    num_rows = num_settings * num_algorithms * num_frequencies * ...
        num_temperatures;
    metrics = local_empty_metric_table(num_rows);
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

        for setting_index = 1:num_settings
            setting = settings(setting_index);
            design = designs{setting_index};
            for algorithm_index = 1:num_algorithms
                algorithm_id = algorithm_ids{algorithm_index};
                for frequency_index = 1:num_frequencies
                    metric_row = metric_row + 1;
                    w = design.filters.(algorithm_id)(:, frequency_index);
                    valid_filter = all(isfinite(w)) && norm(w) > 0;
                    metrics.Setting(metric_row) = string(setting.label);
                    metrics.SettingIndex(metric_row) = setting_index;
                    metrics.Algorithm(metric_row) = ...
                        string(display_names{algorithm_index});
                    metrics.TemperatureC(metric_row) = temperature_c;
                    metrics.FrequencyHz(metric_row) = ...
                        frequencies_hz(frequency_index);
                    metrics.NuB(metric_row) = setting.nu_B(frequency_index);
                    metrics.NuD(metric_row) = setting.nu_D(frequency_index);
                    metrics.ValidFilter(metric_row) = valid_filter;
                    if valid_filter
                        metrics.ControlACdB(metric_row) = local_ac_db(w, ...
                            HB_control(:, :, frequency_index), ...
                            HD_control(:, :, frequency_index));
                        metrics.EvaluationACdB(metric_row) = local_ac_db(w, ...
                            HB_evaluation(:, :, frequency_index), ...
                            HD_evaluation(:, :, frequency_index));
                        if any(strcmp(algorithm_id, {'RPM', 'RACC_PM'}))
                            metrics.ControlNSREdB(metric_row) = ...
                                local_nsre_db(w, ...
                                HB_control(:, :, frequency_index), ...
                                desired_pressure(:, frequency_index));
                        end
                        metrics.EvaluationAEdB(metric_row) = local_ae_db( ...
                            w, HB_evaluation(:, :, frequency_index), ...
                            config.virtual_source_idx, config.target_spl_db);
                    end
                end
            end
        end
        fprintf('Evaluated %.1f deg C (%d/%d).\n', temperature_c, ...
            temperature_index, num_temperatures);
    end

    writetable(metrics, fullfile(output_dir, ...
        'temperature_performance_metrics.csv'));
    summary = local_temperature_summary(metrics, settings, ...
        display_names, frequencies_hz, config.nominal_temperature_c);
    writetable(summary, fullfile(output_dir, ...
        'temperature_performance_summary.csv'));
    figure_paths = local_plot_sweep(summary, settings, display_names, ...
        frequencies_hz, output_dir);

    results = struct();
    results.config = config;
    results.output_dir = output_dir;
    results.active_evaluator = active_evaluator;
    results.active_racc_pm_solver = active_racc_pm_solver;
    results.frequencies_hz = frequencies_hz;
    results.frequency_indices = frequency_indices;
    results.alpha = alpha;
    results.alpha_reference_ac_db = alpha_reference_ac_db;
    results.alpha_reference_filter = alpha_reference_result.w;
    results.settings = settings;
    results.designs = designs;
    results.design_diagnostics = design_diagnostics;
    results.temperature_metrics = metrics;
    results.temperature_summary = summary;
    results.figure_paths = figure_paths;
    save(fullfile(output_dir, 'experiment_results.mat'), 'results', '-v7.3');

    marker_file = fullfile(output_dir, 'RUN_COMPLETE.txt');
    marker_id = fopen(marker_file, 'w');
    if marker_id >= 0
        fprintf(marker_id, ['Completed %s\nSolver: RACC_PM_Sub (no suffix)', ...
            '\nGLS used: no\nrho: %.9g\nmu: %.9g\nalpha margin dB: %.9g\n'], ...
            char(datetime('now')), config.rho, config.mu, ...
            config.alpha_margin_db);
        fclose(marker_id);
    end

    disp(summary);
    fprintf('Completed RACC-PM uncertainty sweep: %s\n', output_dir);
end

function config = local_default_config(response_dir)
    config = struct();
    config.nominal_temperature_c = 22.5;
    config.test_temperatures_c = 20.0:0.1:25.0;
    config.requested_frequencies_hz = [500, 1000, 4000];
    config.relative_uncertainty_levels = ...
        [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.02, 0.05];
    config.include_manuscript_rule = true;
    config.alpha_reference_relative_level = 0.005;
    config.alpha_margin_db = 0;
    config.rho = 10;
    config.mu = 1;
    config.acc_reg_beta = 1e-6;
    config.run_rpm = true;
    config.virtual_source_idx = 13;
    config.target_spl_db = 76;
    config.random_seed = 20260825;
    config.resume = true;
    config.output_dir = '';
    config.output_root = fullfile(response_dir, ...
        'racc_pm_uncertainty_sweep_results');
end

function config = local_apply_overrides(config, overrides)
    if ~isstruct(overrides) || ~isscalar(overrides)
        error('RACCPMRadiusSweep:InvalidConfig', ...
            'user_config must be a scalar structure.');
    end
    names = fieldnames(overrides);
    for k = 1:numel(names)
        if ~isfield(config, names{k})
            error('RACCPMRadiusSweep:UnknownConfigField', ...
                'Unknown configuration field: %s.', names{k});
        end
        config.(names{k}) = overrides.(names{k});
    end
end

function local_validate_config(config)
    validateattributes(config.requested_frequencies_hz, {'numeric'}, ...
        {'vector', 'real', 'finite', 'positive', 'nonempty'});
    validateattributes(config.test_temperatures_c, {'numeric'}, ...
        {'vector', 'real', 'finite', 'nonempty'});
    validateattributes(config.relative_uncertainty_levels, {'numeric'}, ...
        {'vector', 'real', 'finite', 'positive', 'nonempty'});
    validateattributes(config.alpha_reference_relative_level, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(config.alpha_margin_db, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(config.rho, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(config.mu, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(config.acc_reg_beta, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(config.virtual_source_idx, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(config.target_spl_db, {'numeric'}, ...
        {'scalar', 'real', 'finite'});
    if ~islogical(config.include_manuscript_rule) || ...
            ~isscalar(config.include_manuscript_rule)
        error('RACCPMRadiusSweep:InvalidManuscriptFlag', ...
            'include_manuscript_rule must be a scalar logical.');
    end
    if ~islogical(config.run_rpm) || ~isscalar(config.run_rpm)
        error('RACCPMRadiusSweep:InvalidRPMFlag', ...
            'run_rpm must be a scalar logical.');
    end
    if ~islogical(config.resume) || ~isscalar(config.resume)
        error('RACCPMRadiusSweep:InvalidResumeFlag', ...
            'resume must be a scalar logical.');
    end
end

function output_dir = local_prepare_output_directory(config)
    if strlength(string(config.output_dir)) == 0
        timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        output_dir = fullfile(config.output_root, ['run_' timestamp]);
    else
        output_dir = char(string(config.output_dir));
    end
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
end

function indices = local_frequency_indices(all_frequencies_hz, requested_hz)
    indices = zeros(size(requested_hz));
    for k = 1:numel(requested_hz)
        [difference, indices(k)] = min(abs( ...
            all_frequencies_hz - requested_hz(k)));
        if difference > 1e-9
            error('RACCPMRadiusSweep:FrequencyUnavailable', ...
                'Requested frequency %.6g Hz is unavailable.', ...
                requested_hz(k));
        end
    end
end

function [epsilon_B, epsilon_D] = local_relative_radii(HB, HD, nu)
    num_frequencies = size(HB, 3);
    epsilon_B = zeros(1, num_frequencies);
    epsilon_D = zeros(1, num_frequencies);
    for k = 1:num_frequencies
        epsilon_B(k) = nu * norm(HB(:, :, k), 'fro');
        epsilon_D(k) = nu * norm(HD(:, :, k), 'fro');
    end
end

function settings = local_build_settings(config, HB, HD)
    levels = unique(config.relative_uncertainty_levels(:).', 'stable');
    num_constant = numel(levels);
    num_settings = num_constant + double(config.include_manuscript_rule);
    template = struct('id', '', 'label', '', 'type', '', ...
        'nominal_level', NaN, 'epsilon_B', [], 'epsilon_D', [], ...
        'nu_B', [], 'nu_D', []);
    settings = repmat(template, num_settings, 1);
    for k = 1:num_constant
        nu = levels(k);
        [epsilon_B, epsilon_D] = local_relative_radii(HB, HD, nu);
        settings(k).id = ['nu_' local_number_tag(nu)];
        settings(k).label = sprintf('nu=%.4g', nu);
        settings(k).type = 'common-relative';
        settings(k).nominal_level = nu;
        settings(k).epsilon_B = epsilon_B;
        settings(k).epsilon_D = epsilon_D;
        settings(k).nu_B = repmat(nu, 1, size(HB, 3));
        settings(k).nu_D = repmat(nu, 1, size(HD, 3));
    end
    if config.include_manuscript_rule
        index = num_settings;
        num_frequencies = size(HB, 3);
        epsilon_B = zeros(1, num_frequencies);
        epsilon_D = zeros(1, num_frequencies);
        nu_B = zeros(1, num_frequencies);
        nu_D = zeros(1, num_frequencies);
        for k = 1:num_frequencies
            RB = HB(:, :, k)' * HB(:, :, k);
            RD = HD(:, :, k)' * HD(:, :, k);
            epsilon_B(k) = 0.01 * sqrt(norm(RB, 'fro'));
            epsilon_D(k) = 0.01 * sqrt(norm(RD, 'fro'));
            nu_B(k) = epsilon_B(k) / max( ...
                norm(HB(:, :, k), 'fro'), realmin);
            nu_D(k) = epsilon_D(k) / max( ...
                norm(HD(:, :, k), 'fro'), realmin);
        end
        settings(index).id = 'manuscript_rule';
        settings(index).label = 'manuscript rule';
        settings(index).type = 'manuscript-frequency-dependent';
        settings(index).epsilon_B = epsilon_B;
        settings(index).epsilon_D = epsilon_D;
        settings(index).nu_B = nu_B;
        settings(index).nu_D = nu_D;
    end
end

function para = local_matched_para(epsilon_B, epsilon_D, loading_ratio)
    para = struct();
    para.full_racc = struct();
    para.full_racc.eta = struct('B', epsilon_B, 'D', epsilon_D);
    para.full_racc.diagonal_loading_ratio = loading_ratio;
end

function [w_all, diagnostics] = local_design_acc_reg( ...
        HB, HD, beta, frequencies_hz)
    [~, L, num_frequencies] = size(HB);
    w_all = zeros(L, num_frequencies);
    FrequencyHz = frequencies_hz(:);
    Beta = repmat(beta, num_frequencies, 1);
    DenominatorLoading = nan(num_frequencies, 1);
    RawConditionNumber = nan(num_frequencies, 1);
    EffectiveConditionNumber = nan(num_frequencies, 1);
    GeneralizedEigenvalue = nan(num_frequencies, 1);
    RelativeResidual = nan(num_frequencies, 1);
    for k = 1:num_frequencies
        HBi = HB(:, :, k);
        HDi = HD(:, :, k);
        common_scale = max(norm([HBi; HDi], 'fro'), realmin);
        HBs = HBi / common_scale;
        HDs = HDi / common_scale;
        RB = HBs' * HBs;
        RD = HDs' * HDs;
        RB = (RB + RB') / 2;
        RD = (RD + RD') / 2;
        lambda_max_D = max(real(eig(RD)));
        delta = beta * lambda_max_D;
        B = (RD + delta * eye(L));
        B = (B + B') / 2;
        [R, chol_flag] = chol(B);
        if chol_flag ~= 0
            error('RACCPMRadiusSweep:ACCRegNotPositiveDefinite', ...
                'ACC-Reg denominator is not positive definite at %.0f Hz.', ...
                frequencies_hz(k));
        end
        C = R' \ (RB / R);
        C = (C + C') / 2;
        [U, D] = eig(C);
        [lambda, selected_index] = max(real(diag(D)));
        w = R \ U(:, selected_index);
        w = w / max(norm(w), realmin);
        residual = norm(RB * w - lambda * B * w) / max( ...
            norm(RB * w) + abs(lambda) * norm(B * w), realmin);
        w_all(:, k) = w;
        DenominatorLoading(k) = delta * common_scale^2;
        RawConditionNumber(k) = cond(RD);
        EffectiveConditionNumber(k) = cond(B);
        GeneralizedEigenvalue(k) = lambda;
        RelativeResidual(k) = residual;
    end
    diagnostics = table(FrequencyHz, Beta, DenominatorLoading, ...
        RawConditionNumber, EffectiveConditionNumber, ...
        GeneralizedEigenvalue, RelativeResidual);
end

function status = local_design_status(design, algorithm_id, ...
        frequency_index, valid_filter)
    if strcmp(algorithm_id, 'RACC_PM')
        if isempty(design.racc_pm_result)
            status = "Exception";
        else
            status = string(design.racc_pm_result.status(frequency_index));
        end
    elseif strcmp(algorithm_id, 'RPM') && ...
            strlength(string(design.rpm_error)) > 0
        status = "Exception";
    elseif valid_filter
        status = "Solved";
    else
        status = "Invalid filter";
    end
end

function table_out = local_empty_design_table(num_rows)
    table_out = table(strings(num_rows, 1), strings(num_rows, 1), ...
        nan(num_rows, 1), nan(num_rows, 1), nan(num_rows, 1), ...
        nan(num_rows, 1), nan(num_rows, 1), nan(num_rows, 1), ...
        false(num_rows, 1), nan(num_rows, 1), strings(num_rows, 1), ...
        nan(num_rows, 1), nan(num_rows, 1), nan(num_rows, 1), ...
        'VariableNames', {'Setting', 'Algorithm', 'FrequencyHz', 'NuB', ...
        'NuD', 'EpsilonB', 'EpsilonD', 'Alpha', 'ValidFilter', ...
        'FilterNorm', 'Status', 'ControlACdB', 'ControlNSREdB', ...
        'ControlAEdB'});
end

function table_out = local_empty_metric_table(num_rows)
    table_out = table(strings(num_rows, 1), nan(num_rows, 1), ...
        strings(num_rows, 1), nan(num_rows, 1), nan(num_rows, 1), ...
        nan(num_rows, 1), nan(num_rows, 1), false(num_rows, 1), ...
        nan(num_rows, 1), nan(num_rows, 1), nan(num_rows, 1), ...
        nan(num_rows, 1), ...
        'VariableNames', {'Setting', 'SettingIndex', 'Algorithm', ...
        'TemperatureC', 'FrequencyHz', 'NuB', 'NuD', 'ValidFilter', ...
        'ControlACdB', 'EvaluationACdB', 'ControlNSREdB', ...
        'EvaluationAEdB'});
end

function summary = local_temperature_summary(metrics, settings, ...
        display_names, frequencies_hz, nominal_temperature_c)
    num_rows = numel(settings) * numel(display_names) * ...
        numel(frequencies_hz);
    Setting = strings(num_rows, 1);
    SettingIndex = nan(num_rows, 1);
    Algorithm = strings(num_rows, 1);
    FrequencyHz = nan(num_rows, 1);
    NuB = nan(num_rows, 1);
    NuD = nan(num_rows, 1);
    NominalEvaluationACdB = nan(num_rows, 1);
    MeanEvaluationACdB = nan(num_rows, 1);
    MinimumEvaluationACdB = nan(num_rows, 1);
    EvaluationACRangedB = nan(num_rows, 1);
    NominalControlNSREdB = nan(num_rows, 1);
    MeanControlNSREdB = nan(num_rows, 1);
    MaximumControlNSREdB = nan(num_rows, 1);
    NominalEvaluationAEdB = nan(num_rows, 1);
    MeanEvaluationAEdB = nan(num_rows, 1);
    MaximumEvaluationAEdB = nan(num_rows, 1);
    row = 0;
    for setting_index = 1:numel(settings)
        for algorithm_index = 1:numel(display_names)
            for frequency_index = 1:numel(frequencies_hz)
                row = row + 1;
                selected = metrics.SettingIndex == setting_index & ...
                    metrics.Algorithm == string(display_names{algorithm_index}) & ...
                    metrics.FrequencyHz == frequencies_hz(frequency_index);
                temperatures = metrics.TemperatureC(selected);
                ac = metrics.EvaluationACdB(selected);
                nsre = metrics.ControlNSREdB(selected);
                ae = metrics.EvaluationAEdB(selected);
                nominal = abs(temperatures - nominal_temperature_c) < 1e-9;
                Setting(row) = string(settings(setting_index).label);
                SettingIndex(row) = setting_index;
                Algorithm(row) = string(display_names{algorithm_index});
                FrequencyHz(row) = frequencies_hz(frequency_index);
                NuB(row) = settings(setting_index).nu_B(frequency_index);
                NuD(row) = settings(setting_index).nu_D(frequency_index);
                NominalEvaluationACdB(row) = local_nominal(ac, nominal);
                MeanEvaluationACdB(row) = mean(ac, 'omitnan');
                MinimumEvaluationACdB(row) = min(ac, [], 'omitnan');
                EvaluationACRangedB(row) = ...
                    max(ac, [], 'omitnan') - min(ac, [], 'omitnan');
                NominalControlNSREdB(row) = local_nominal(nsre, nominal);
                MeanControlNSREdB(row) = mean(nsre, 'omitnan');
                MaximumControlNSREdB(row) = max(nsre, [], 'omitnan');
                NominalEvaluationAEdB(row) = local_nominal(ae, nominal);
                MeanEvaluationAEdB(row) = mean(ae, 'omitnan');
                MaximumEvaluationAEdB(row) = max(ae, [], 'omitnan');
            end
        end
    end
    summary = table(Setting, SettingIndex, Algorithm, FrequencyHz, NuB, ...
        NuD, NominalEvaluationACdB, MeanEvaluationACdB, ...
        MinimumEvaluationACdB, EvaluationACRangedB, ...
        NominalControlNSREdB, MeanControlNSREdB, ...
        MaximumControlNSREdB, NominalEvaluationAEdB, ...
        MeanEvaluationAEdB, MaximumEvaluationAEdB);
end

function value = local_nominal(values, nominal_mask)
    selected = values(nominal_mask);
    selected = selected(isfinite(selected));
    if isempty(selected)
        value = NaN;
    else
        value = selected(1);
    end
end

function paths = local_plot_sweep(summary, settings, display_names, ...
        frequencies_hz, output_dir)
    metric_specs = { ...
        'NominalEvaluationACdB', 'Nominal evaluation AC (dB)', ...
            'nominal_ac_vs_uncertainty'; ...
        'MinimumEvaluationACdB', 'Minimum evaluation AC (dB)', ...
            'minimum_ac_vs_uncertainty'; ...
        'NominalControlNSREdB', 'Nominal control NSRE (dB)', ...
            'nominal_nsre_vs_uncertainty'; ...
        'MaximumControlNSREdB', 'Maximum control NSRE (dB)', ...
            'maximum_nsre_vs_uncertainty'; ...
        'NominalEvaluationAEdB', 'Nominal evaluation AE (dB)', ...
            'nominal_ae_vs_uncertainty'};
    colors = lines(numel(display_names));
    line_styles = {'-', '--', '-.', ':'};
    markers = {'o', 's', '^', 'd'};
    paths = struct();
    constant_setting_indices = find(strcmp({settings.type}, ...
        'common-relative'));
    manuscript_setting_index = find(strcmp({settings.type}, ...
        'manuscript-frequency-dependent'), 1);

    for metric_index = 1:size(metric_specs, 1)
        variable_name = metric_specs{metric_index, 1};
        y_label = metric_specs{metric_index, 2};
        base_name = metric_specs{metric_index, 3};
        figure_handle = figure('Visible', 'off', 'Color', 'w', ...
            'Position', [50, 80, 1380, 430]);
        layout = tiledlayout(figure_handle, 1, numel(frequencies_hz), ...
            'TileSpacing', 'compact', 'Padding', 'compact');
        legend_handles = gobjects(1, numel(display_names));

        for frequency_index = 1:numel(frequencies_hz)
            ax = nexttile(layout);
            hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
            for algorithm_index = 1:numel(display_names)
                selected = ismember(summary.SettingIndex, ...
                    constant_setting_indices) & ...
                    summary.Algorithm == string(display_names{algorithm_index}) & ...
                    summary.FrequencyHz == frequencies_hz(frequency_index);
                x = summary.NuB(selected);
                y = summary.(variable_name)(selected);
                [x, order] = sort(x);
                y = y(order);
                h = semilogx(ax, x, y, ...
                    'LineStyle', line_styles{algorithm_index}, ...
                    'Marker', markers{algorithm_index}, ...
                    'Color', colors(algorithm_index, :), ...
                    'LineWidth', 1.6, 'MarkerSize', 5, ...
                    'DisplayName', display_names{algorithm_index});
                if frequency_index == 1
                    legend_handles(algorithm_index) = h;
                end
                if ~isempty(manuscript_setting_index)
                    paper = summary.SettingIndex == ...
                        manuscript_setting_index & ...
                        summary.Algorithm == ...
                        string(display_names{algorithm_index}) & ...
                        summary.FrequencyHz == frequencies_hz(frequency_index);
                    semilogx(ax, summary.NuB(paper), ...
                        summary.(variable_name)(paper), ...
                        'LineStyle', 'none', 'Marker', 'p', ...
                        'MarkerSize', 9, 'LineWidth', 1.4, ...
                        'Color', colors(algorithm_index, :), ...
                        'HandleVisibility', 'off');
                end
            end
            xlabel(ax, 'Relative ATF radius \nu_B');
            if frequency_index == 1
                ylabel(ax, y_label);
            end
            title(ax, sprintf('%d Hz', frequencies_hz(frequency_index)), ...
                'FontWeight', 'normal');
            set(ax, 'FontName', 'Times New Roman', 'FontSize', 11, ...
                'LineWidth', 1);
        end
        legend_object = legend(legend_handles, display_names, ...
            'Orientation', 'horizontal', ...
            'NumColumns', numel(display_names), ...
            'FontName', 'Times New Roman', 'FontSize', 10);
        legend_object.Layout.Tile = 'south';
        drawnow;
        png_path = fullfile(output_dir, [base_name '.png']);
        pdf_path = fullfile(output_dir, [base_name '.pdf']);
        fig_path = fullfile(output_dir, [base_name '.fig']);
        exportgraphics(figure_handle, png_path, 'Resolution', 300);
        exportgraphics(figure_handle, pdf_path, 'ContentType', 'vector');
        savefig(figure_handle, fig_path);
        close(figure_handle);
        paths.(base_name) = struct('png', png_path, ...
            'pdf', pdf_path, 'fig', fig_path);
    end
end

function ac_db = local_ac_db(w, HB, HD)
    bright_energy = norm(HB * w)^2;
    dark_energy = norm(HD * w)^2;
    MB = size(HB, 1);
    MD = size(HD, 1);
    contrast = (MD * bright_energy) / max(MB * dark_energy, realmin);
    ac_db = 10 * log10(max(real(contrast), realmin));
end

function nsre_db = local_nsre_db(w, HB, desired_pressure)
    error_energy = norm(HB * w - desired_pressure)^2;
    target_energy = norm(desired_pressure)^2;
    nsre_db = 10 * log10(max(error_energy / max( ...
        target_energy, realmin), realmin));
end

function ae_db = local_ae_db(w, HB, virtual_source_index, target_spl_db)
    if virtual_source_index > size(HB, 2)
        error('RACCPMRadiusSweep:VirtualSourceOutOfRange', ...
            'virtual_source_idx exceeds the loudspeaker count.');
    end
    reference = zeros(size(w));
    reference(virtual_source_index) = 1;
    target_pressure = 20e-6 * 10^(target_spl_db / 20);
    current_rms = sqrt(mean(abs(HB * w).^2));
    calibrated_w = w * target_pressure / max(current_rms, realmin);
    RB = HB' * HB;
    energy_ratio = real(reference' * RB * reference) / max( ...
        real(calibrated_w' * RB * calibrated_w), realmin);
    effort_ratio = real(calibrated_w' * calibrated_w);
    ae_db = 10 * log10(max(effort_ratio * energy_ratio, realmin));
end

function tag = local_number_tag(value)
    tag = strrep(sprintf('%.9g', value), '.', 'p');
    tag = strrep(tag, '-', 'm');
    tag = strrep(tag, '+', 'p');
end
