function outputs = run_sicer_vast_temperature_experiment(output_dir)
%RUN_SICER_VAST_TEMPERATURE_EXPERIMENT Temperature-informed VAST benchmark.
%   This standalone experiment preserves the authoritative ten-algorithm
%   result and adds two adaptive SICER-VAST endpoint benchmarks:
%       SICER_VAST_V1  - DFT-domain VAST rank 1 (ACC endpoint)
%       SICER_VAST_VL  - DFT-domain VAST rank L (PM endpoint)
%
%   At every target temperature, the nominal 22.5 C control RIRs are SICER
%   corrected, VAST is redesigned from the corrected control model, and the
%   filters are evaluated on the true independent evaluation ATFs at that
%   temperature. Ground-truth redesigned VAST endpoints are calculated only
%   as a reproduction diagnostic and are not added as competing algorithms.

    response_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(response_dir);
    if ~strcmpi(pwd, project_root)
        error('SICERVAST:WrongWorkingDirectory', ...
            'Run this function from the project root: %s.', project_root);
    end
    addpath(response_dir);
    addpath(genpath(fullfile(project_root, 'src')));
    addpath(genpath(fullfile(project_root, 'data')));
    addpath(fullfile(project_root, 'src', 'evaluations'), '-begin');

    expected_evaluator = fullfile(project_root, 'src', 'evaluations', ...
        'evaluate_performance.m');
    resolved_evaluator = which('evaluate_performance');
    if ~strcmpi(resolved_evaluator, expected_evaluator)
        error('SICERVAST:WrongEvaluator', ...
            'Resolved evaluator is %s; expected %s.', ...
            resolved_evaluator, expected_evaluator);
    end

    base_result_file = fullfile(response_dir, ...
        'racc_pm_parameter_comparison_results', ...
        'rho10_nu0p01_all_run_20260827_002614', ...
        'performance_matrices_20260827_021428.mat');
    base_analysis_dir = fullfile(fileparts(base_result_file), 'analysis');
    monitor_nsre_file = fullfile(base_analysis_dir, ...
        'monitor_nsre_postprocess.mat');
    if ~isfile(base_result_file) || ~isfile(monitor_nsre_file)
        error('SICERVAST:MissingAuthoritativeSource', ...
            'The authoritative result or monitor NSRE sidecar is missing.');
    end

    base = load(base_result_file, 'all_filters', 'algorithms_to_test', ...
        'algorithm_wall_seconds', 'val_results', 'para', 'freq_params');
    monitor = load(monitor_nsre_file, 'results');
    if nargin < 1 || isempty(output_dir)
        run_id = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        output_dir = fullfile(response_dir, ...
            'sicer_vast_temperature_results', ['run_' run_id]);
    else
        output_dir = char(string(output_dir));
        run_id = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    end
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    figure_dir = fullfile(output_dir, 'figures');
    snapshot_dir = fullfile(output_dir, 'perceptual_filter_snapshots');
    if ~exist(figure_dir, 'dir'), mkdir(figure_dir); end
    if ~exist(snapshot_dir, 'dir'), mkdir(snapshot_dir); end

    diary_file = fullfile(output_dir, 'matlab_diary.log');
    diary('off');
    diary(diary_file);
    diary_cleanup = onCleanup(@() diary('off'));
    start_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
    marker = fopen(fullfile(output_dir, 'RUN_STARTED.txt'), 'w');
    if marker >= 0
        fprintf(marker, 'Started: %s\n', char(start_time));
        fprintf(marker, 'MATLAB PID: %d\n', feature('getpid'));
        fprintf(marker, 'Base result: %s\n', base_result_file);
        fclose(marker);
    end
    fprintf('SICER-VAST temperature experiment started: %s\n', ...
        char(start_time));
    fprintf('Output: %s\n', output_dir);
    fprintf('Evaluator: %s\n', resolved_evaluator);
    fprintf('Authoritative base: %s\n', base_result_file);

    data_dir = fullfile(project_root, 'data', 'SimulateRIR', 'temperature');
    nominal_file = get_data_filename(data_dir, 'temperature', 22.5);
    nominal = load(nominal_file, 'IR_BZ', 'IR_DZ', 'ATF_BZ', 'ATF_DZ');
    desired_file = fullfile(project_root, 'data', 'ATF_desired_plane.mat');
    desired_loaded = load(desired_file, 'ATF_desired_plane');
    desired_control = desired_loaded.ATF_desired_plane;
    freq_params = base.para.freq_params;
    frequencies = freq_params.target_freqs(:).';
    temperatures = base.para.temperature_vector_celsius(:);
    num_temperatures = numel(temperatures);
    num_frequencies = numel(frequencies);
    loudspeakers = size(nominal.ATF_BZ.ctrl, 2);
    bright_points = size(nominal.IR_BZ.ctrl, 1);

    geometry = load(fullfile(project_root, 'data', 'arrayGeometry', ...
        'array_layout.mat'), 'roomArray');
    desired_monitor = local_plane_wave(geometry.roomArray, ...
        geometry.roomArray.BZ_eval, freq_params);

    endpoint_algorithms = {'SICER_VAST_V1', 'SICER_VAST_VL'};
    endpoint_fields = {'V1', sprintf('V%d', loudspeakers)};
    gt_algorithms = {'GT_VAST_V1', 'GT_VAST_VL'};
    metric_names = {'AC', 'RawMonitorNSRE', 'NSRE', ...
        'ContrastCalibratedError', 'AE', 'Planarity', 'BZ_SPL_Std'};
    val_results = local_initialize_results(endpoint_algorithms, ...
        metric_names, num_temperatures, num_frequencies);
    gt_results = local_initialize_results(gt_algorithms, metric_names, ...
        num_temperatures, num_frequencies);
    filters_by_temperature = struct();
    gt_filters_by_temperature = struct();
    for a = 1:numel(endpoint_algorithms)
        filters_by_temperature.(endpoint_algorithms{a}) = complex(zeros( ...
            loudspeakers, num_frequencies, num_temperatures));
        gt_filters_by_temperature.(gt_algorithms{a}) = complex(zeros( ...
            loudspeakers, num_frequencies, num_temperatures));
    end

    base_algorithms = base.algorithms_to_test;
    base_bz_spl_std = struct();
    for a = 1:numel(base_algorithms)
        base_bz_spl_std.(base_algorithms{a}) = ...
            nan(num_temperatures, num_frequencies);
    end

    corrected_control_atfs = struct();
    selected_snapshot_temperatures = [20, 22.5, 25];
    for k = 1:numel(selected_snapshot_temperatures)
        corrected_control_atfs.(local_temperature_field( ...
            selected_snapshot_temperatures(k))) = [];
    end

    correction_wall_seconds = nan(num_temperatures, 1);
    correction_cpu_seconds = nan(num_temperatures, 1);
    vast_wall_seconds = nan(num_temperatures, 1);
    vast_cpu_seconds = nan(num_temperatures, 1);
    gt_vast_wall_seconds = nan(num_temperatures, 1);
    beta_values = nan(num_temperatures, 1);
    correction_nmse_db = nan(num_temperatures, 1);
    correction_magnitude_mae_db = nan(num_temperatures, 1);
    min_dark_rcond = nan(num_temperatures, 1);
    max_full_rank_residual = nan(num_temperatures, 1);

    combined_nominal_rir = cat(1, nominal.IR_BZ.ctrl, nominal.IR_DZ.ctrl);
    old_sound_speed = temp2speed(22.5);
    vast_para = struct();
    vast_para.vast.ranks = [1, loudspeakers];
    vast_para.vast.mu = 1;
    vast_para.vast.dark_loading_ratio = 0;
    vast_para.vast.fail_rcond_threshold = 0;

    total_timer = tic;
    for temperature_index = 1:num_temperatures
        temperature = temperatures(temperature_index);
        fprintf('\n[%d/%d] Temperature %.1f C\n', temperature_index, ...
            num_temperatures, temperature);
        operating_file = get_data_filename(data_dir, 'temperature', ...
            temperature);
        operating = load(operating_file, 'IR_BZ', 'IR_DZ', ...
            'ATF_BZ', 'ATF_DZ');
        new_sound_speed = temp2speed(temperature);

        correction_cpu_start = cputime;
        correction_timer = tic;
        [corrected_combined, correction_info] = sicer_correct_rir( ...
            combined_nominal_rir, old_sound_speed, new_sound_speed, ...
            freq_params.fs, struct('lowpass_order', 100, ...
            'apply_antialias', true, 'chunk_size', 512));
        correction_wall_seconds(temperature_index) = toc(correction_timer);
        correction_cpu_seconds(temperature_index) = ...
            cputime - correction_cpu_start;
        beta_values(temperature_index) = correction_info.beta;

        corrected_BZ_rir = corrected_combined(1:bright_points, :, :);
        corrected_DZ_rir = corrected_combined(bright_points+1:end, :, :);
        corrected_BZ_atf = compute_atf(corrected_BZ_rir, freq_params);
        corrected_DZ_atf = compute_atf(corrected_DZ_rir, freq_params);

        [correction_nmse_db(temperature_index), ...
            correction_magnitude_mae_db(temperature_index)] = ...
            local_atf_correction_error(corrected_BZ_atf, corrected_DZ_atf, ...
            operating.ATF_BZ.ctrl, operating.ATF_DZ.ctrl, frequencies);

        snapshot_index = find(abs(selected_snapshot_temperatures - ...
            temperature) < 1e-10, 1);
        if ~isempty(snapshot_index)
            snapshot_field = local_temperature_field(temperature);
            corrected_control_atfs.(snapshot_field) = struct( ...
                'BZ', corrected_BZ_atf, 'DZ', corrected_DZ_atf);
        end

        vast_cpu_start = cputime;
        vast_timer = tic;
        [sicer_filters, vast_info] = vast_dft_endpoints( ...
            corrected_BZ_atf, corrected_DZ_atf, desired_control, vast_para);
        vast_wall_seconds(temperature_index) = toc(vast_timer);
        vast_cpu_seconds(temperature_index) = cputime - vast_cpu_start;
        min_dark_rcond(temperature_index) = min(vast_info.dark_rcond);
        max_full_rank_residual(temperature_index) = ...
            max(vast_info.full_rank_relative_residual);

        gt_timer = tic;
        [gt_filters, ~] = vast_dft_endpoints(operating.ATF_BZ.ctrl, ...
            operating.ATF_DZ.ctrl, desired_control, vast_para);
        gt_vast_wall_seconds(temperature_index) = toc(gt_timer);

        steering = precompute_steering_matrix( ...
            geometry.roomArray.BZ_eval, frequencies, new_sound_speed);
        for endpoint_index = 1:numel(endpoint_algorithms)
            algorithm = endpoint_algorithms{endpoint_index};
            endpoint_field = endpoint_fields{endpoint_index};
            weights = sicer_filters.(endpoint_field);
            filters_by_temperature.(algorithm)(:, :, temperature_index) = ...
                weights;
            metrics = local_evaluate(weights, operating, corrected_BZ_atf, ...
                desired_control, desired_monitor, steering, ...
                strcmp(algorithm, 'SICER_VAST_VL'));
            val_results = local_assign_metrics(val_results, algorithm, ...
                temperature_index, metrics, metric_names);

            gt_algorithm = gt_algorithms{endpoint_index};
            gt_weights = gt_filters.(endpoint_field);
            gt_filters_by_temperature.(gt_algorithm)(:, :, ...
                temperature_index) = gt_weights;
            gt_metrics = local_evaluate(gt_weights, operating, ...
                operating.ATF_BZ.ctrl, desired_control, desired_monitor, ...
                steering, strcmp(gt_algorithm, 'GT_VAST_VL'));
            gt_results = local_assign_metrics(gt_results, gt_algorithm, ...
                temperature_index, gt_metrics, metric_names);
        end

        for base_index = 1:numel(base_algorithms)
            algorithm = base_algorithms{base_index};
            base_bz_spl_std.(algorithm)(temperature_index, :) = ...
                local_bz_spl_std(base.all_filters.(algorithm), ...
                operating.ATF_BZ.eval);
        end
        clear corrected_combined corrected_BZ_rir corrected_DZ_rir;
        fprintf(['  beta=%.8f, correction %.2f s, VAST %.2f s, ', ...
            'median ATF NMSE %.2f dB\n'], beta_values(temperature_index), ...
            correction_wall_seconds(temperature_index), ...
            vast_wall_seconds(temperature_index), ...
            correction_nmse_db(temperature_index));
    end
    total_wall_seconds = toc(total_timer);

    combined_algorithms = [base_algorithms, endpoint_algorithms];
    combined_results = base.val_results;
    for a = 1:numel(base_algorithms)
        algorithm = base_algorithms{a};
        if local_is_target_referenced(algorithm)
            combined_results.(algorithm).NSRE = ...
                monitor.results.raw_monitor_nsre_db.(algorithm);
            combined_results.(algorithm).ContrastCalibratedError = ...
                nan(size(combined_results.(algorithm).NSRE));
        else
            combined_results.(algorithm).NSRE = ...
                nan(size(monitor.results.raw_monitor_nsre_db.(algorithm)));
            combined_results.(algorithm).ContrastCalibratedError = ...
                monitor.results.aligned_monitor_nsre_db.(algorithm);
        end
        combined_results.(algorithm).BZ_SPL_Std = ...
            base_bz_spl_std.(algorithm);
    end
    for a = 1:numel(endpoint_algorithms)
        algorithm = endpoint_algorithms{a};
        combined_results.(algorithm).AC = val_results.(algorithm).AC;
        combined_results.(algorithm).NSRE = val_results.(algorithm).NSRE;
        combined_results.(algorithm).ContrastCalibratedError = ...
            val_results.(algorithm).ContrastCalibratedError;
        combined_results.(algorithm).AE = val_results.(algorithm).AE;
        combined_results.(algorithm).Planarity = ...
            val_results.(algorithm).Planarity;
        combined_results.(algorithm).BZ_SPL_Std = ...
            val_results.(algorithm).BZ_SPL_Std;
    end

    runtime = table(temperatures, beta_values, correction_wall_seconds, ...
        correction_cpu_seconds, vast_wall_seconds, vast_cpu_seconds, ...
        gt_vast_wall_seconds, correction_nmse_db, ...
        correction_magnitude_mae_db, min_dark_rcond, ...
        max_full_rank_residual, 'VariableNames', ...
        {'TemperatureC', 'Beta', 'SICERWallSeconds', 'SICERCPUSeconds', ...
         'VASTWallSeconds', 'VASTCPUSeconds', 'GTVASTWallSeconds', ...
         'MedianATFComplexNMSEdB', 'MedianATFMagnitudeMAEdB', ...
         'MinDarkRcond', 'MaxFullRankResidual'});
    writetable(runtime, fullfile(output_dir, ...
        'sicer_vast_runtime_and_correction.csv'));

    summary = local_combined_summary(combined_results, ...
        combined_algorithms, temperatures, frequencies, ...
        base.algorithm_wall_seconds, correction_wall_seconds, ...
        vast_wall_seconds);
    writetable(summary, fullfile(output_dir, ...
        'combined_100_4000Hz_temperature_summary.csv'));
    gt_gap_summary = local_gt_gap_summary(val_results, gt_results, ...
        endpoint_algorithms, gt_algorithms, metric_names, frequencies);
    writetable(gt_gap_summary, fullfile(output_dir, ...
        'sicer_vs_ground_truth_gap_summary.csv'));

    figure_paths = local_plot_combined(combined_results, ...
        combined_algorithms, temperatures, frequencies, figure_dir);
    snapshot_files = local_save_perceptual_snapshots(snapshot_dir, ...
        selected_snapshot_temperatures, temperatures, ...
        endpoint_algorithms, filters_by_temperature, val_results, ...
        corrected_control_atfs, base.para, base_result_file);

    result_file = fullfile(output_dir, 'sicer_vast_temperature_results.mat');
    experiment = struct();
    experiment.run_id = run_id;
    experiment.output_dir = output_dir;
    experiment.base_result_file = base_result_file;
    experiment.evaluator = resolved_evaluator;
    experiment.temperatures = temperatures;
    experiment.frequencies = frequencies;
    experiment.endpoint_algorithms = endpoint_algorithms;
    experiment.gt_algorithms = gt_algorithms;
    experiment.vast_parameters = vast_para.vast;
    experiment.sicer_parameters = struct('nominal_temperature_celsius', ...
        22.5, 'lowpass_order', 100, 'antialias', true, ...
        'chunk_size', 512);
    experiment.val_results = val_results;
    experiment.gt_results = gt_results;
    experiment.filters_by_temperature = filters_by_temperature;
    experiment.gt_filters_by_temperature = gt_filters_by_temperature;
    experiment.combined_algorithms = combined_algorithms;
    experiment.combined_results = combined_results;
    experiment.runtime = runtime;
    experiment.summary = summary;
    experiment.gt_gap_summary = gt_gap_summary;
    experiment.figure_paths = figure_paths;
    experiment.snapshot_files = snapshot_files;
    experiment.total_wall_seconds = total_wall_seconds;
    save(result_file, 'experiment', '-v7.3');

    end_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
    marker = fopen(fullfile(output_dir, 'RUN_COMPLETE.txt'), 'w');
    if marker >= 0
        fprintf(marker, 'Completed: %s\n', char(end_time));
        fprintf(marker, 'Result: %s\n', result_file);
        fprintf(marker, 'Wall seconds: %.6f\n', total_wall_seconds);
        fclose(marker);
    end
    fprintf('\nSICER-VAST experiment completed: %s\n', char(end_time));
    fprintf('Total wall time: %.2f s\n', total_wall_seconds);
    fprintf('Result: %s\n', result_file);
    disp(summary);

    outputs = experiment;
    outputs.result_file = result_file;
    clear diary_cleanup;
end

function results = local_initialize_results(algorithms, metrics, rows, columns)
    results = struct();
    for a = 1:numel(algorithms)
        for m = 1:numel(metrics)
            results.(algorithms{a}).(metrics{m}) = nan(rows, columns);
        end
    end
end

function results = local_assign_metrics(results, algorithm, row, metrics, names)
    for m = 1:numel(names)
        results.(algorithm).(names{m})(row, :) = metrics.(names{m});
    end
end

function metrics = local_evaluate(weights, operating, design_HB, ...
        desired_control, desired_monitor, steering, is_target_referenced)
    num_frequencies = size(weights, 2);
    metrics = struct('AC', nan(1, num_frequencies), ...
        'RawMonitorNSRE', nan(1, num_frequencies), ...
        'NSRE', nan(1, num_frequencies), ...
        'ContrastCalibratedError', nan(1, num_frequencies), ...
        'AE', nan(1, num_frequencies), ...
        'Planarity', nan(1, num_frequencies), ...
        'BZ_SPL_Std', nan(1, num_frequencies));
    for f = 1:num_frequencies
        q = weights(:, f);
        HB_eval = operating.ATF_BZ.eval(:, :, f);
        HD_eval = operating.ATF_DZ.eval(:, :, f);
        metrics.AC(f) = real(calculate_AC(q, HB_eval, HD_eval));
        metrics.AE(f) = real(calculate_AE(q, HB_eval, 13, 76));
        metrics.Planarity(f) = calculate_planarity(q, HB_eval, ...
            steering(:, :, f));
        metrics.RawMonitorNSRE(f) = calculate_NSRE(q, HB_eval, ...
            desired_monitor(:, f));

        design_pressure = design_HB(:, :, f) * q;
        desired = desired_control(:, f);
        cross_term = design_pressure' * desired;
        if abs(cross_term) <= realmin
            phase_factor = 1;
        else
            phase_factor = exp(1i * angle(cross_term));
        end
        gain = (norm(desired) / max(norm(design_pressure), realmin)) * ...
            phase_factor;
        if is_target_referenced
            metrics.NSRE(f) = metrics.RawMonitorNSRE(f);
        else
            metrics.ContrastCalibratedError(f) = calculate_NSRE( ...
                gain * q, HB_eval, desired_monitor(:, f));
        end
    end
    metrics.BZ_SPL_Std = local_bz_spl_std(weights, ...
        operating.ATF_BZ.eval);
end

function spatial_std = local_bz_spl_std(weights, HB_eval)
    num_frequencies = size(weights, 2);
    spatial_std = nan(1, num_frequencies);
    for f = 1:num_frequencies
        pressure = HB_eval(:, :, f) * weights(:, f);
        levels = 20 * log10(max(abs(pressure), realmin) / 20e-6);
        spatial_std(f) = std(levels, 1);
    end
end

function desired = local_plane_wave(room_array, positions, freq_params)
    source_position = room_array.s(13, :);
    propagation = [2.6, 2.0, 1.6] - source_position;
    propagation = propagation / norm(propagation);
    wave_numbers = (2 * pi .* freq_params.target_freqs(:)) / 343;
    phase_delay = positions * (wave_numbers .* propagation).';
    desired = exp(-1i * phase_delay);
end

function [median_nmse_db, median_magnitude_mae_db] = ...
        local_atf_correction_error(corrected_BZ, corrected_DZ, ...
        true_BZ, true_DZ, frequencies)
    corrected = cat(1, corrected_BZ, corrected_DZ);
    truth = cat(1, true_BZ, true_DZ);
    paths = size(corrected, 1) * size(corrected, 2);
    corrected = reshape(corrected, paths, size(corrected, 3));
    truth = reshape(truth, paths, size(truth, 3));
    sample_indices = round(linspace(1, paths, min(96, paths)));
    band = frequencies >= 100 & frequencies <= 4000;
    error_value = corrected(sample_indices, band) - ...
        truth(sample_indices, band);
    nmse = 10 * log10(sum(abs(error_value).^2, 2) ./ ...
        max(sum(abs(truth(sample_indices, band)).^2, 2), realmin));
    magnitude_mae = mean(abs(20 * log10(max(abs( ...
        corrected(sample_indices, band)), realmin)) - ...
        20 * log10(max(abs(truth(sample_indices, band)), realmin))), 2);
    median_nmse_db = median(nmse);
    median_magnitude_mae_db = median(magnitude_mae);
end

function summary = local_combined_summary(results, algorithms, ...
        temperatures, frequencies, base_times, correction_times, vast_times)
    band = frequencies >= 100 & frequencies <= 4000;
    nominal_index = local_temperature_index(temperatures, 22.5);
    rows = numel(algorithms);
    Algorithm = string(algorithms(:));
    UpdateMode = strings(rows, 1);
    KeyParameters = strings(rows, 1);
    NominalACdB = nan(rows, 1); SweepMeanACdB = nan(rows, 1);
    WorstACdB = nan(rows, 1); NominalNSREdB = nan(rows, 1);
    SweepMeanNSREdB = nan(rows, 1); WorstNSREdB = nan(rows, 1);
    NominalPlanarityPct = nan(rows, 1); SweepMeanPlanarityPct = nan(rows, 1);
    WorstPlanarityPct = nan(rows, 1); NominalBZSPLStd_dB = nan(rows, 1);
    SweepMeanBZSPLStd_dB = nan(rows, 1); WorstBZSPLStd_dB = nan(rows, 1);
    NominalAEdB = nan(rows, 1); SweepMeanAEdB = nan(rows, 1);
    WorstAEdB = nan(rows, 1); RuntimeSeconds = nan(rows, 1);
    RuntimeBasis = strings(rows, 1);

    for a = 1:rows
        algorithm = algorithms{a};
        ac = mean(results.(algorithm).AC(:, band), 2, 'omitnan');
        nsre = mean(results.(algorithm).NSRE(:, band), 2, 'omitnan');
        planarity = mean(results.(algorithm).Planarity(:, band), 2, 'omitnan');
        bz_std = mean(results.(algorithm).BZ_SPL_Std(:, band), 2, 'omitnan');
        ae = mean(results.(algorithm).AE(:, band), 2, 'omitnan');
        NominalACdB(a) = ac(nominal_index);
        SweepMeanACdB(a) = mean(ac, 'omitnan');
        WorstACdB(a) = min(ac, [], 'omitnan');
        NominalNSREdB(a) = nsre(nominal_index);
        SweepMeanNSREdB(a) = mean(nsre, 'omitnan');
        WorstNSREdB(a) = max(nsre, [], 'omitnan');
        NominalPlanarityPct(a) = planarity(nominal_index);
        SweepMeanPlanarityPct(a) = mean(planarity, 'omitnan');
        WorstPlanarityPct(a) = min(planarity, [], 'omitnan');
        NominalBZSPLStd_dB(a) = bz_std(nominal_index);
        SweepMeanBZSPLStd_dB(a) = mean(bz_std, 'omitnan');
        WorstBZSPLStd_dB(a) = max(bz_std, [], 'omitnan');
        NominalAEdB(a) = ae(nominal_index);
        SweepMeanAEdB(a) = mean(ae, 'omitnan');
        WorstAEdB(a) = max(ae, [], 'omitnan');

        if startsWith(algorithm, 'SICER_')
            UpdateMode(a) = "Temperature-informed redesign";
            if endsWith(algorithm, 'V1')
                KeyParameters(a) = "V=1, mu=1";
            else
                KeyParameters(a) = "V=L, mu=1";
            end
            RuntimeSeconds(a) = mean(correction_times + vast_times);
            RuntimeBasis(a) = "SICER + VAST per update";
        else
            UpdateMode(a) = "Fixed filter at 22.5 C";
            KeyParameters(a) = "Authoritative rho=10, nu=0.01 run";
            if isfield(base_times, algorithm)
                RuntimeSeconds(a) = base_times.(algorithm);
            end
            RuntimeBasis(a) = "One 159-bin design";
        end
    end
    summary = table(Algorithm, UpdateMode, KeyParameters, NominalACdB, ...
        SweepMeanACdB, WorstACdB, NominalNSREdB, SweepMeanNSREdB, ...
        WorstNSREdB, NominalPlanarityPct, SweepMeanPlanarityPct, ...
        WorstPlanarityPct, NominalBZSPLStd_dB, ...
        SweepMeanBZSPLStd_dB, WorstBZSPLStd_dB, NominalAEdB, ...
        SweepMeanAEdB, WorstAEdB, RuntimeSeconds, RuntimeBasis);
end

function summary = local_gt_gap_summary(sicer, gt, algorithms, ...
        gt_algorithms, metrics, frequencies)
    band = frequencies >= 100 & frequencies <= 4000;
    rows = numel(algorithms) * numel(metrics);
    Algorithm = strings(rows, 1); Metric = strings(rows, 1);
    MeanAbsoluteGap = nan(rows, 1); MaxAbsoluteGap = nan(rows, 1);
    row = 0;
    for a = 1:numel(algorithms)
        for m = 1:numel(metrics)
            row = row + 1;
            values = mean(sicer.(algorithms{a}).(metrics{m})(:, band), ...
                2, 'omitnan');
            reference = mean(gt.(gt_algorithms{a}).(metrics{m})(:, band), ...
                2, 'omitnan');
            gap = values - reference;
            Algorithm(row) = string(algorithms{a});
            Metric(row) = string(metrics{m});
            MeanAbsoluteGap(row) = mean(abs(gap), 'omitnan');
            MaxAbsoluteGap(row) = max(abs(gap), [], 'omitnan');
        end
    end
    summary = table(Algorithm, Metric, MeanAbsoluteGap, MaxAbsoluteGap);
end

function paths = local_plot_combined(results, algorithms, temperatures, ...
        frequencies, output_dir)
    band = frequencies >= 100 & frequencies <= 4000;
    plot_metrics = {'AC', 'NSRE', 'Planarity', 'BZ_SPL_Std'};
    ylabels = {'AC (dB)', 'Standard monitor NSRE (dB)', ...
        'Planarity (%)', 'BZ SPL spatial std. (dB)'};
    colors = lines(numel(algorithms));
    fig = figure('Color', 'w', 'Position', [100, 100, 1240, 780]);
    layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
        'Padding', 'compact');
    handles = gobjects(numel(algorithms), 1);
    for panel = 1:numel(plot_metrics)
        ax = nexttile(layout);
        hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        for a = 1:numel(algorithms)
            values = mean(results.(algorithms{a}).(plot_metrics{panel})(:, band), ...
                2, 'omitnan');
            line_width = 1.0;
            line_style = '-';
            if startsWith(algorithms{a}, 'SICER_')
                line_width = 2.2;
                line_style = '--';
            end
            handles(a) = plot(ax, temperatures, values, ...
                'Color', colors(a, :), 'LineWidth', line_width, ...
                'LineStyle', line_style);
        end
        xline(ax, 22.5, ':k', 'HandleVisibility', 'off');
        xlabel(ax, 'Temperature (^{\circ}C)');
        ylabel(ax, ylabels{panel});
        set(ax, 'FontName', 'Times New Roman', 'FontSize', 10);
    end
    display_names = cellfun(@local_display_name, algorithms, ...
        'UniformOutput', false);
    legend_handle = legend(handles, display_names, 'NumColumns', 4, ...
        'Orientation', 'horizontal', 'Location', 'southoutside');
    if isprop(legend_handle, 'Layout')
        legend_handle.Layout.Tile = 'south';
    end
    base = fullfile(output_dir, ...
        'combined_sicer_vast_temperature_metrics');
    paths = {[base '.png'], [base '.pdf'], [base '.fig']};
    exportgraphics(fig, paths{1}, 'Resolution', 300);
    exportgraphics(fig, paths{2}, 'ContentType', 'vector');
    savefig(fig, paths{3});
    close(fig);
end

function tf = local_is_target_referenced(algorithm)
    tf = ismember(string(algorithm), ["PM", "ACC_PM", "WCRPM", ...
        "RPM", "RACC_PM", "RACC_PM_Sub", "RACC_PM_Subpro", ...
        "RACC_PM_GLS", "SICER_VAST_VL", "GT_VAST_VL"]);
end

function snapshot_files = local_save_perceptual_snapshots(output_dir, ...
        requested_temperatures, temperatures, algorithms, filters, ...
        val_results, corrected_atfs, para, base_result_file)
    snapshot_files = cell(numel(requested_temperatures), 1);
    for k = 1:numel(requested_temperatures)
        temperature = requested_temperatures(k);
        temperature_index = local_temperature_index(temperatures, temperature);
        all_filters = struct();
        snapshot_val_results = struct();
        for a = 1:numel(algorithms)
            algorithm = algorithms{a};
            all_filters.(algorithm) = filters.(algorithm)(:, :, ...
                temperature_index);
            snapshot_val_results.(algorithm) = val_results.(algorithm);
        end
        calibration_field = local_temperature_field(temperature);
        perceptual_config = struct();
        perceptual_config.algorithm_fields = string(algorithms);
        perceptual_config.display_algorithms = string(algorithms);
        perceptual_config.listening_algorithms = string(algorithms);
        perceptual_config.display_names = [ ...
            "SICER-VAST (V=1)", "SICER-VAST (V=L)"];
        perceptual_config.design_temperature_celsius = temperature;
        perceptual_config.calibration_temperature_celsius = temperature;
        perceptual_config.method_note = ...
            "Temperature-informed SICER correction and VAST redesign";
        snapshot_files{k} = fullfile(output_dir, sprintf( ...
            'sicer_vast_filters_T%05.2fC.mat', temperature));
        snapshot = struct();
        snapshot.all_filters = all_filters;
        snapshot.algorithms_to_test = algorithms;
        snapshot.para = para;
        snapshot.val_results = snapshot_val_results;
        snapshot.calibration_ATF_BZ_ctrl = ...
            corrected_atfs.(calibration_field).BZ;
        snapshot.calibration_temperature_celsius = temperature;
        snapshot.filter_temperature_celsius = temperature;
        snapshot.perceptual_config = perceptual_config;
        snapshot.source_result_file = base_result_file;
        save(snapshot_files{k}, '-struct', 'snapshot', '-v7.3');
    end
end

function index = local_temperature_index(temperatures, requested)
    [distance, index] = min(abs(temperatures - requested));
    if distance > 1e-10
        error('SICERVAST:TemperatureMissing', ...
            'Temperature %.2f C is unavailable.', requested);
    end
end

function field = local_temperature_field(temperature)
    field = matlab.lang.makeValidName(sprintf('T_%05.2fC', temperature));
end

function name = local_display_name(identifier)
    switch identifier
        case 'ACC_Reg', name = 'ACC-Reg';
        case 'ACC_PM', name = 'ACC-PM';
        case 'wcRACC', name = 'WCRACC';
        case 'NoCT_WCRACC', name = 'NoCT-WCRACC';
        case 'Full_WCRACC', name = 'Full-WCRACC';
        case 'POTDC_RACC', name = 'POTDC-RACC';
        case 'RACC_PM_Subpro', name = 'RACC-PM';
        case 'SICER_VAST_V1', name = 'SICER-VAST (V=1)';
        case 'SICER_VAST_VL', name = 'SICER-VAST (V=L)';
        otherwise, name = strrep(identifier, '_', '-');
    end
end
