function output = run_full_cross_term_racc_experiment(run_mode, run_dir, parallel_workers)
%RUN_FULL_CROSS_TERM_RACC_EXPERIMENT Validate and run the Full-RACC ablation.
%
% Examples (run from the project root):
%   run_full_cross_term_racc_experiment('validate')
%   run_full_cross_term_racc_experiment('smoke-real')
%   run_full_cross_term_racc_experiment('frequency')
%   run_full_cross_term_racc_experiment('full')
%   run_full_cross_term_racc_experiment('full', [], 4) % explicit workers
%
% Modes:
%   validate   - small complex system; raw/decomposed SDP equivalence test.
%   smoke-real - one real 1-kHz bin with a coarse bisection tolerance.
%   frequency  - real data at 0.5, 1, 2, and 4 kHz.
%   full       - all 159 bins from 100 Hz to 8 kHz.
%
% The nominal filter is designed at 22.5 deg C. Absolute ATF uncertainty
% radii cover every simulated temperature from 22.2 to 22.8 deg C:
%
%   eta_Z(f) = max_T ||H_Z(f,T) - H_Z(f,22.5)||_F.
%
% Therefore the conventional ablation uses exactly gamma_Z(f)=eta_Z(f)^2.
% Each real-data Full-RACC frequency is checkpointed independently. To
% resume an interrupted run, pass the existing run directory as run_dir.
% The decomposed LMI is intentionally locked to MOSEK; using SDPT3/SeDuMi
% at the 96-control-point, 48-source scale can run out of memory.

    if nargin < 1 || isempty(run_mode)
        run_mode = 'validate';
    end
    if nargin < 3
        parallel_workers = [];
    end
    run_mode = lower(char(string(run_mode)));

    response_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(response_dir);
    original_dir = pwd;
    cleanup = onCleanup(@() cd(original_dir)); %#ok<NASGU>
    cd(project_root);
    addpath(genpath(fullfile(project_root, 'src')));
    addpath(response_dir);

    % AGENTS.md requires checking the duplicate function before any MATLAB
    % experiment. The evaluations implementation must win path resolution.
    resolved_evaluation = which('evaluate_performance');
    expected_evaluation = fullfile(project_root, 'src', 'evaluations', ...
        'evaluate_performance.m');
    if ~strcmpi(resolved_evaluation, expected_evaluation)
        error('FullRACC:WrongEvaluatePerformance', ...
            ['MATLAB resolves evaluate_performance to:\n  %s\nExpected:\n  %s\n', ...
             'Fix the path order before running the experiment.'], ...
            resolved_evaluation, expected_evaluation);
    end
    fprintf('evaluate_performance resolves to: %s\n', resolved_evaluation);

    if nargin < 2 || isempty(run_dir)
        timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        run_dir = fullfile(response_dir, 'full_cross_term_racc_results', ...
            [strrep(run_mode, '-', '_') '_' timestamp]);
    elseif ~isfolder(run_dir)
        error('FullRACC:ResumeDirectoryMissing', ...
            'The requested resume directory does not exist: %s', run_dir);
    end
    if ~isfolder(run_dir)
        mkdir(run_dir);
    end

    switch run_mode
        case 'validate'
            output = local_validate_formulations(run_dir);
        case {'smoke-real', 'frequency', 'full'}
            output = local_run_temperature_ablation( ...
                run_mode, run_dir, project_root, parallel_workers);
        otherwise
            error('FullRACC:UnknownRunMode', ...
                'Unknown mode "%s". Use validate, smoke-real, frequency, or full.', ...
                run_mode);
    end

    output.run_mode = run_mode;
    output.run_dir = run_dir;
    save(fullfile(run_dir, 'experiment_output.mat'), 'output', '-v7.3');
    fprintf('Full-RACC experiment output: %s\n', run_dir);
end

function output = local_validate_formulations(run_dir)
    rng(20260820, 'twister');
    MB = 4;
    MD = 5;
    L = 3;
    HB = (randn(MB, L) + 1i * randn(MB, L)) / sqrt(2);
    HD = (randn(MD, L) + 1i * randn(MD, L)) / sqrt(2);
    etaB = 0.02 * norm(HB, 'fro');
    etaD = 0.02 * norm(HD, 'fro');

    para = struct();
    para.full_racc.eta = struct('B', etaB, 'D', etaD);
    para.full_racc.solver = 'mosek';
    para.full_racc.bisection_tolerance_db = 0.10;
    para.full_racc.max_bisection_iterations = 14;
    para.full_racc.randomization_trials = 200;
    para.full_racc.random_seed = 20260820;
    para.full_racc.mosek_threads = 4;

    para.full_racc.solution_method = 'bisection';
    para.full_racc.formulation = 'raw';
    raw = FullCrossTermRACC(HB, HD, para);
    para.full_racc.formulation = 'decomposed';
    decomposed = FullCrossTermRACC(HB, HD, para);
    para.full_racc.solution_method = 'normalized';
    normalized = FullCrossTermRACC(HB, HD, para);

    relaxation_difference_db = abs( ...
        raw.relaxation_ac_db(1) - decomposed.relaxation_ac_db(1));
    tolerance_db = para.full_racc.bisection_tolerance_db + 1e-6;
    if relaxation_difference_db > tolerance_db
        error('FullRACC:DecompositionValidationFailed', ...
            ['Raw and decomposed relaxation levels differ by %.4f dB, ', ...
             'which exceeds %.4f dB.'], relaxation_difference_db, tolerance_db);
    end
    normalized_difference_db = abs( ...
        normalized.relaxation_ac_db(1) - decomposed.relaxation_ac_db(1));
    if normalized_difference_db > tolerance_db
        error('FullRACC:NormalizationValidationFailed', ...
            ['Single-SDP normalization and bisection differ by %.4f dB, ', ...
             'which exceeds %.4f dB.'], normalized_difference_db, tolerance_db);
    end

    % Construct perturbations that attain the closed-form rank-one bounds.
    w = decomposed.w(:, 1);
    [formula_db, formula_alpha, formula_details] = ...
        full_racc_worst_case_contrast(w, HB, HD, etaB, etaD);
    [attained_alpha, attained_details] = ...
        local_construct_attaining_perturbations(w, HB, HD, etaB, etaD);
    closed_form_error = abs(formula_alpha - attained_alpha);
    if closed_form_error > 1e-8 * max(1, abs(formula_alpha))
        error('FullRACC:WorstCaseFormulaValidationFailed', ...
            'The constructed perturbation does not attain the closed-form bound.');
    end

    summary = table( ...
        string({'Raw'; 'Decomposed'; 'Decomposed'}), ...
        string({'Bisection'; 'Bisection'; 'Normalized'}), ...
        [raw.relaxation_ac_db(1); decomposed.relaxation_ac_db(1); ...
        normalized.relaxation_ac_db(1)], ...
        [raw.recovered_ac_db(1); decomposed.recovered_ac_db(1); ...
        normalized.recovered_ac_db(1)], ...
        [raw.rank_fraction(1); decomposed.rank_fraction(1); ...
        normalized.rank_fraction(1)], ...
        [raw.solve_time_seconds(1); decomposed.solve_time_seconds(1); ...
        normalized.solve_time_seconds(1)], ...
        string({raw.status{1}; decomposed.status{1}; normalized.status{1}}), ...
        'VariableNames', {'Formulation', 'SolutionMethod', 'RelaxationACdB', ...
        'RecoveredWorstCaseACdB', 'RankFraction', 'SolveTimeSeconds', 'Status'});
    writetable(summary, fullfile(run_dir, 'raw_vs_decomposed_validation.csv'));

    output = struct();
    output.raw = raw;
    output.decomposed = decomposed;
    output.normalized = normalized;
    output.summary = summary;
    output.relaxation_difference_db = relaxation_difference_db;
    output.normalized_difference_db = normalized_difference_db;
    output.closed_form = struct('ac_db', formula_db, ...
        'alpha', formula_alpha, 'details', formula_details, ...
        'attained_alpha', attained_alpha, ...
        'attained_details', attained_details, ...
        'absolute_error', closed_form_error);

    fprintf('Raw/decomposed relaxation difference: %.6g dB\n', ...
        relaxation_difference_db);
    fprintf('Normalized/bisection relaxation difference: %.6g dB\n', ...
        normalized_difference_db);
    fprintf('Recovered exact ATF-ball worst-case AC: %.3f dB\n', formula_db);
    disp(summary);
end

function output = local_run_temperature_ablation( ...
        run_mode, run_dir, project_root, requested_parallel_workers)
    data_dir = fullfile(project_root, 'data', 'SimulateRIR', 'temperature');
    para_data = load(fullfile(data_dir, 'para.mat'), 'para');
    freq_params = para_data.para.freq_params;
    all_frequencies = freq_params.target_freqs(:).';

    switch run_mode
        case 'smoke-real'
            requested_frequencies = 1000;
            tolerance_db = 1.0;
            max_iterations = 7;
            randomization_trials = 50;
        case 'frequency'
            requested_frequencies = [500, 1000, 2000, 4000];
            tolerance_db = 0.20;
            max_iterations = 12;
            randomization_trials = 300;
        case 'full'
            requested_frequencies = all_frequencies;
            tolerance_db = 0.20;
            max_iterations = 12;
            randomization_trials = 300;
    end

    frequency_indices = local_nearest_indices(all_frequencies, requested_frequencies);
    frequencies = all_frequencies(frequency_indices);
    [parallel_workers, mosek_threads] = local_parallel_configuration( ...
        requested_parallel_workers, numel(frequencies), run_mode);
    fprintf('Full-RACC runtime: %d frequency worker(s), %d MOSEK thread(s) per worker.\n', ...
        parallel_workers, mosek_threads);
    nominal_temperature = 22.5;
    uncertainty_temperatures = 22.2:0.1:22.8;

    nominal_file = get_data_filename(data_dir, 'temperature', nominal_temperature);
    nominal_data = load(nominal_file, 'ATF_BZ', 'ATF_DZ');
    HB_nominal = nominal_data.ATF_BZ.ctrl(:, :, frequency_indices);
    HD_nominal = nominal_data.ATF_DZ.ctrl(:, :, frequency_indices);
    clear nominal_data

    [etaB, etaD] = local_temperature_bounds( ...
        data_dir, nominal_temperature, uncertainty_temperatures, ...
        frequency_indices, HB_nominal, HD_nominal);

    experiment_config = struct( ...
        'nominal_temperature_celsius', nominal_temperature, ...
        'uncertainty_temperatures_celsius', uncertainty_temperatures, ...
        'frequency_indices', frequency_indices, ...
        'frequencies_hz', frequencies, ...
        'etaB', etaB, 'etaD', etaD, ...
        'gammaB', etaB.^2, 'gammaD', etaD.^2, ...
        'solver', 'mosek', 'formulation', 'decomposed', ...
        'solution_method', 'normalized', ...
        'bisection_tolerance_db', tolerance_db, ...
        'max_bisection_iterations', max_iterations, ...
        'randomization_trials', randomization_trials);
    config_file = fullfile(run_dir, 'experiment_config.mat');
    if isfile(config_file)
        previous = load(config_file, 'experiment_config');
        if ~isequaln(previous.experiment_config, experiment_config)
            error('FullRACC:ResumeConfigurationMismatch', ...
                'The requested run does not match the saved checkpoint configuration.');
        end
    else
        save(config_file, 'experiment_config');
    end

    acc_w = ACC(HB_nominal, HD_nominal);
    matched_para.full_racc.eta = struct('B', etaB, 'D', etaD);
    conventional = ConventionalRACCMatched( ...
        HB_nominal, HD_nominal, matched_para);

    num_frequencies = numel(frequencies);
    L = size(HB_nominal, 2);
    full_results = cell(1, num_frequencies);
    pending_indices = zeros(1, num_frequencies);
    pending_count = 0;
    for k = 1:num_frequencies
        checkpoint_file = local_checkpoint_file(run_dir, k, frequencies(k));
        if isfile(checkpoint_file)
            checkpoint = load(checkpoint_file, 'full_result');
            full_results{k} = checkpoint.full_result;
            fprintf('Loaded checkpoint for %.1f Hz.\n', frequencies(k));
        else
            pending_count = pending_count + 1;
            pending_indices(pending_count) = k;
        end
    end
    pending_indices = pending_indices(1:pending_count);

    if parallel_workers > 1 && pending_count > 1
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool('local', parallel_workers);
        end
        active_workers = min([parallel_workers, pool.NumWorkers, pending_count]);
        pending_results = cell(1, pending_count);
        parfor (job = 1:pending_count, active_workers)
            k = pending_indices(job);
            pending_results{job} = local_solve_frequency_checkpoint( ...
                k, frequencies(k), run_dir, HB_nominal(:, :, k), ...
                HD_nominal(:, :, k), etaB(k), etaD(k), tolerance_db, ...
                max_iterations, randomization_trials, mosek_threads);
        end
        for job = 1:pending_count
            full_results{pending_indices(job)} = pending_results{job};
        end
    else
        for job = 1:pending_count
            k = pending_indices(job);
            full_results{k} = local_solve_frequency_checkpoint( ...
                k, frequencies(k), run_dir, HB_nominal(:, :, k), ...
                HD_nominal(:, :, k), etaB(k), etaD(k), tolerance_db, ...
                max_iterations, randomization_trials, mosek_threads);
        end
    end

    full_w = zeros(L, num_frequencies);
    full_relaxation_db = nan(1, num_frequencies);
    full_recovered_db = nan(1, num_frequencies);
    full_rank_fraction = nan(1, num_frequencies);
    full_solve_time = nan(1, num_frequencies);
    full_status = repmat({''}, 1, num_frequencies);
    full_recovery_source = repmat({''}, 1, num_frequencies);

    for k = 1:num_frequencies
        full_result = full_results{k};

        full_w(:, k) = full_result.w(:, 1);
        full_relaxation_db(k) = full_result.relaxation_ac_db(1);
        full_recovered_db(k) = full_result.recovered_ac_db(1);
        full_rank_fraction(k) = full_result.rank_fraction(1);
        full_solve_time(k) = full_result.solve_time_seconds(1);
        full_status{k} = full_result.status{1};
        if isfield(full_result, 'recovery_source')
            full_recovery_source{k} = full_result.recovery_source{1};
        end
    end

    filters = struct('ACC', acc_w, ...
        'ConventionalRACC', conventional.w, 'FullRACC', full_w);
    desired_data = load(fullfile(project_root, 'data', ...
        'ATF_desired_plane.mat'), 'ATF_desired_plane');
    desired_pressure = desired_data.ATF_desired_plane(:, frequency_indices);
    evaluation_table = local_evaluate_temperature_metrics( ...
        data_dir, uncertainty_temperatures, frequency_indices, ...
        frequencies, filters, desired_pressure);
    writetable(evaluation_table, fullfile(run_dir, ...
        'temperature_metrics_across_algorithms.csv'));
    variation_table = local_temperature_variation_summary( ...
        evaluation_table, frequencies, fieldnames(filters));
    writetable(variation_table, fullfile(run_dir, ...
        'temperature_metric_variation_summary.csv'));

    acc_worst_db = local_score_filters(acc_w, HB_nominal, HD_nominal, etaB, etaD);
    conventional_worst_db = local_score_filters( ...
        conventional.w, HB_nominal, HD_nominal, etaB, etaD);
    full_worst_db = local_score_filters( ...
        full_w, HB_nominal, HD_nominal, etaB, etaD);
    diagnostic_table = table(frequencies(:), etaB(:), etaD(:), ...
        etaB(:).^2, etaD(:).^2, acc_worst_db(:), ...
        conventional_worst_db(:), full_worst_db(:), ...
        full_relaxation_db(:), full_recovered_db(:), ...
        full_rank_fraction(:), full_solve_time(:), string(full_status(:)), ...
        string(full_recovery_source(:)), ...
        'VariableNames', {'FrequencyHz', 'EtaB', 'EtaD', 'GammaB', 'GammaD', ...
        'ACCExactWorstCaseACdB', 'ConventionalExactWorstCaseACdB', ...
        'FullExactWorstCaseACdB', 'FullSDPRelaxationACdB', ...
        'FullRecoveredACdB', 'FullRankFraction', ...
        'FullSolveTimeSeconds', 'FullStatus', 'FullRecoverySource'});
    writetable(diagnostic_table, fullfile(run_dir, ...
        'full_racc_diagnostics.csv'));

    output = struct();
    output.config = experiment_config;
    output.filters = filters;
    output.conventional = conventional;
    output.evaluation_table = evaluation_table;
    output.variation_table = variation_table;
    output.diagnostic_table = diagnostic_table;
    output.runtime = struct('parallel_workers', parallel_workers, ...
        'mosek_threads_per_worker', mosek_threads, ...
        'estimated_memory_gb_per_worker', 10);
    disp(diagnostic_table);
end

function full_result = local_solve_frequency_checkpoint( ...
        k, frequency, run_dir, HB, HD, etaB, etaD, tolerance_db, ...
        max_iterations, randomization_trials, mosek_threads)
    full_para = struct();
    full_para.full_racc.eta = struct('B', etaB, 'D', etaD);
    full_para.full_racc.formulation = 'decomposed';
    full_para.full_racc.solution_method = 'normalized';
    full_para.full_racc.solver = 'mosek';
    full_para.full_racc.mosek_threads = mosek_threads;
    full_para.full_racc.bisection_tolerance_db = tolerance_db;
    full_para.full_racc.max_bisection_iterations = max_iterations;
    full_para.full_racc.randomization_trials = randomization_trials;
    full_para.full_racc.random_seed = 20260820 + k;
    full_para.full_racc.refinement_iterations = 300;
    full_result = FullCrossTermRACC(HB, HD, full_para);
    checkpoint_file = local_checkpoint_file(run_dir, k, frequency);
    save(checkpoint_file, 'full_result', '-v7.3');
end

function checkpoint_file = local_checkpoint_file(run_dir, k, frequency)
    checkpoint_file = fullfile(run_dir, sprintf( ...
        'frequency_%03d_%07.1fHz.mat', k, frequency));
end

function [workers, mosek_threads] = local_parallel_configuration( ...
        requested_workers, num_frequencies, run_mode)
    physical_cores = feature('numcores');
    if strcmp(run_mode, 'smoke-real')
        workers = 1;
    elseif ~isempty(requested_workers)
        validateattributes(requested_workers, {'numeric'}, ...
            {'scalar', 'integer', 'positive', 'finite'});
        workers = requested_workers;
    else
        % A measured 96x48 decomposed solve uses roughly 8--9 GB on the
        % current workstation. Reserve 8 GB for the client/OS and budget
        % 10 GB per solver worker. The default cap avoids oversubscribing
        % MOSEK even on machines with many logical cores.
        try
            [~, system_memory] = memory;
            available_bytes = system_memory.PhysicalMemory.Available;
            memory_workers = floor(max(available_bytes - 8 * 1024^3, 0) / ...
                (10 * 1024^3));
        catch
            memory_workers = 1;
        end
        cpu_workers = max(1, floor(physical_cores / 2));
        workers = max(1, min([4, memory_workers, cpu_workers]));
    end
    workers = max(1, min(workers, num_frequencies));
    mosek_threads = max(1, min(4, floor(physical_cores / workers)));
end

function indices = local_nearest_indices(all_frequencies, requested)
    indices = zeros(size(requested));
    for k = 1:numel(requested)
        [~, indices(k)] = min(abs(all_frequencies - requested(k)));
    end
    indices = unique(indices, 'stable');
end

function [etaB, etaD] = local_temperature_bounds(data_dir, nominal_temperature, ...
        temperatures, frequency_indices, HB_nominal, HD_nominal)
    num_frequencies = numel(frequency_indices);
    etaB = zeros(1, num_frequencies);
    etaD = zeros(1, num_frequencies);
    for temperature = temperatures
        if abs(temperature - nominal_temperature) < 1e-10
            continue;
        end
        data_file = get_data_filename(data_dir, 'temperature', temperature);
        current = load(data_file, 'ATF_BZ', 'ATF_DZ');
        HB = current.ATF_BZ.ctrl(:, :, frequency_indices);
        HD = current.ATF_DZ.ctrl(:, :, frequency_indices);
        for k = 1:num_frequencies
            etaB(k) = max(etaB(k), norm(HB(:, :, k) - HB_nominal(:, :, k), 'fro'));
            etaD(k) = max(etaD(k), norm(HD(:, :, k) - HD_nominal(:, :, k), 'fro'));
        end
        clear current HB HD
    end
end

function table_out = local_evaluate_temperature_metrics(data_dir, temperatures, ...
        frequency_indices, frequencies, filters, desired_pressure)
    algorithm_names = fieldnames(filters);
    num_rows = numel(temperatures) * numel(frequencies) * numel(algorithm_names);
    temperature_column = zeros(num_rows, 1);
    frequency_column = zeros(num_rows, 1);
    algorithm_column = strings(num_rows, 1);
    ac_column = zeros(num_rows, 1);
    nsre_column = zeros(num_rows, 1);
    ae_column = zeros(num_rows, 1);
    row = 0;

    for temperature = temperatures
        data_file = get_data_filename(data_dir, 'temperature', temperature);
        current = load(data_file, 'ATF_BZ', 'ATF_DZ');
        HB_eval = current.ATF_BZ.eval(:, :, frequency_indices);
        HD_eval = current.ATF_DZ.eval(:, :, frequency_indices);
        HB_ctrl = current.ATF_BZ.ctrl(:, :, frequency_indices);
        for a = 1:numel(algorithm_names)
            algorithm = algorithm_names{a};
            w = filters.(algorithm);
            for k = 1:numel(frequencies)
                row = row + 1;
                temperature_column(row) = temperature;
                frequency_column(row) = frequencies(k);
                algorithm_column(row) = string(algorithm);
                ac_column(row) = real(calculate_AC( ...
                    w(:, k), HB_eval(:, :, k), HD_eval(:, :, k)));
                nsre_column(row) = real(calculate_NSRE( ...
                    w(:, k), HB_ctrl(:, :, k), desired_pressure(:, k)));
                ae_column(row) = real(calculate_AE( ...
                    w(:, k), HB_eval(:, :, k), 13, 76));
            end
        end
        clear current HB_eval HD_eval HB_ctrl
    end

    table_out = table(temperature_column, frequency_column, ...
        algorithm_column, ac_column, nsre_column, ae_column, ...
        'VariableNames', {'TemperatureC', 'FrequencyHz', 'Algorithm', ...
        'ACdB', 'NSREdB', 'AEdB'});
end

function summary = local_temperature_variation_summary( ...
        evaluation_table, frequencies, algorithm_names)
    num_rows = numel(frequencies) * numel(algorithm_names);
    frequency_column = zeros(num_rows, 1);
    algorithm_column = strings(num_rows, 1);
    ac_mean = zeros(num_rows, 1);
    ac_min = zeros(num_rows, 1);
    ac_range = zeros(num_rows, 1);
    nsre_mean = zeros(num_rows, 1);
    nsre_range = zeros(num_rows, 1);
    ae_mean = zeros(num_rows, 1);
    ae_range = zeros(num_rows, 1);
    row = 0;
    for a = 1:numel(algorithm_names)
        algorithm = string(algorithm_names{a});
        for k = 1:numel(frequencies)
            row = row + 1;
            selected = evaluation_table.Algorithm == algorithm & ...
                evaluation_table.FrequencyHz == frequencies(k);
            ac = evaluation_table.ACdB(selected);
            nsre = evaluation_table.NSREdB(selected);
            ae = evaluation_table.AEdB(selected);
            frequency_column(row) = frequencies(k);
            algorithm_column(row) = algorithm;
            ac_mean(row) = mean(ac);
            ac_min(row) = min(ac);
            ac_range(row) = max(ac) - min(ac);
            nsre_mean(row) = mean(nsre);
            nsre_range(row) = max(nsre) - min(nsre);
            ae_mean(row) = mean(ae);
            ae_range(row) = max(ae) - min(ae);
        end
    end
    summary = table(frequency_column, algorithm_column, ac_mean, ac_min, ...
        ac_range, nsre_mean, nsre_range, ae_mean, ae_range, ...
        'VariableNames', {'FrequencyHz', 'Algorithm', 'MeanACdB', ...
        'MinimumACdB', 'ACRangeAcrossTemperaturedB', 'MeanNSREdB', ...
        'NSRERangeAcrossTemperaturedB', 'MeanAEdB', ...
        'AERangeAcrossTemperaturedB'});
end

function scores_db = local_score_filters(w_all, HB, HD, etaB, etaD)
    num_frequencies = size(w_all, 2);
    scores_db = nan(1, num_frequencies);
    for k = 1:num_frequencies
        scores_db(k) = full_racc_worst_case_contrast( ...
            w_all(:, k), HB(:, :, k), HD(:, :, k), etaB(k), etaD(k));
    end
end

function [alpha, details] = local_construct_attaining_perturbations( ...
        w, HB, HD, etaB, etaD)
    w_norm_sq = real(w' * w);
    yB = HB * w;
    yD = HD * w;

    if norm(yB) <= etaB * sqrt(w_norm_sq)
        delta_yB = -yB;
    elseif norm(yB) > 0
        delta_yB = -etaB * sqrt(w_norm_sq) * yB / norm(yB);
    else
        delta_yB = zeros(size(yB));
    end
    if norm(yD) > 0
        delta_yD = etaD * sqrt(w_norm_sq) * yD / norm(yD);
    else
        delta_yD = zeros(size(yD));
        if etaD > 0
            delta_yD(1) = etaD * sqrt(w_norm_sq);
        end
    end

    deltaHB = delta_yB * w' / w_norm_sq;
    deltaHD = delta_yD * w' / w_norm_sq;
    bright_energy = norm((HB + deltaHB) * w)^2;
    dark_energy = norm((HD + deltaHD) * w)^2;
    MB = size(HB, 1);
    MD = size(HD, 1);
    alpha = (MD * bright_energy) / (MB * dark_energy);
    details = struct('deltaHB_fro', norm(deltaHB, 'fro'), ...
        'deltaHD_fro', norm(deltaHD, 'fro'), ...
        'bright_energy', bright_energy, 'dark_energy', dark_energy);
end
