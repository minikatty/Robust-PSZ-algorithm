%% FINAL FAIR ALL-ALGORITHM EVALUATION (NO V2 RACC-PM)
% The RACC-PM uncertainty radius is frozen before this script is run. Every
% algorithm uses the same nominal ATFs, frequency grid, loudspeaker layout,
% physical test ATFs, and matched uncertainty level. Evaluation points and
% off-nominal temperatures are not used to tune alpha or uncertainty.

clear; clc; close all;

%% ========================================================================
%  1. CONFIGURATION & SETUP
%  ========================================================================
project_root = fileparts(fileparts(mfilename('fullpath')));
if ~strcmpi(pwd, project_root)
    error('FairEvaluation:WrongWorkingDirectory', ...
        'Run this script from the project root: %s.', project_root);
end
addpath(fullfile(project_root, 'RQ_response'));
addpath(genpath(fullfile(project_root, 'src')));
addpath(genpath(fullfile(project_root, 'data')));
addpath(fullfile(project_root, 'src', 'evaluations'), '-begin');
% --- Experiment Configuration ---
experiment_modes = {'temperature'};
algorithm_profile = strtrim(getenv('FAIR_ALGORITHM_PROFILE'));
if isempty(algorithm_profile), algorithm_profile = 'original'; end
switch lower(algorithm_profile)
    case 'original'
        algorithms_to_test = {'ACC_Reg', 'PM', 'ACC_PM', 'wcRACC', ...
            'NoCT_WCRACC', 'Full_WCRACC', 'POTDC_RACC', 'RPM', ...
            'RACC_PM_Subpro'};
    case 'racc_pm_with_acc'
        algorithms_to_test = {'ACC', 'RACC_PM_Subpro'};
    case 'all_with_acc'
        algorithms_to_test = {'ACC', 'ACC_Reg', 'PM', 'ACC_PM', ...
            'wcRACC', 'NoCT_WCRACC', 'Full_WCRACC', 'POTDC_RACC', ...
            'RPM', 'RACC_PM_Subpro'};
    otherwise
        error('FairEvaluation:UnknownAlgorithmProfile', ...
            'Unknown FAIR_ALGORITHM_PROFILE: %s.', algorithm_profile);
end
metrics_to_evaluate = {'AC', 'NSRE', 'AE', 'Planarity'};
ctrl_design_env = [22.5, 0.1, 20]; % the enviroment for each case

% Common uncertainty configuration for every revised robust method.
relative_uncertainty_level = local_environment_scalar('FAIR_NU', 5e-3);
racc_pm_rho = local_environment_scalar('FAIR_RHO', 10);
if relative_uncertainty_level <= 0 || racc_pm_rho <= 0
    error('FairEvaluation:InvalidParameterOverride', ...
        'FAIR_NU and FAIR_RHO must both be positive.');
end
diagonal_loading_ratio = 1e-6;
full_wcracc_parallel_workers = 2;
mosek_threads_per_worker = 2;
racc_pm_parallel_workers = 1;
racc_pm_mosek_threads_per_worker = [];
reuse_full_checkpoint_dir = fullfile(project_root, 'results', ...
    'temperature_AnalysisData', ...
    'full_wcracc_checkpoints_20260824_000937');

%% ========================================================================
%  2. MAIN EXECUTION BLOCK
%  ========================================================================
total_timer = tic;
% --- Loop through each experiment mode ---
for mode_idx = 1:length(experiment_modes)
    current_mode = experiment_modes{mode_idx};
    mode_timer = tic;

    fprintf('====== Starting Experiment Mode: %s ======\n', upper(current_mode));
    run_timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    requested_run_dir = getenv('FAIR_RUN_DIR');
    if isempty(requested_run_dir)
        results_dir = fullfile(project_root, 'RQ_response', ...
            'all_algorithms_fair_results', ['run_' run_timestamp]);
    else
        results_dir = requested_run_dir;
        fprintf('Resuming fair-run directory: %s\n', results_dir);
    end
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    diary(fullfile(results_dir, 'matlab_diary.log'));
    start_marker_id = fopen(fullfile(results_dir, 'RUN_STARTED.txt'), 'w');
    if start_marker_id >= 0
        fprintf(start_marker_id, 'Started: %s\n', ...
            char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
        fprintf(start_marker_id, 'MATLAB PID: %d\n', feature('getpid'));
        fprintf(start_marker_id, 'Script: %s\n', mfilename('fullpath'));
        fclose(start_marker_id);
    end
    design_checkpoint_dir = fullfile(results_dir, 'design_checkpoints');
    if ~exist(design_checkpoint_dir, 'dir'), mkdir(design_checkpoint_dir); end

    requested_checkpoint_dir = getenv('FULL_WCRACC_CHECKPOINT_DIR');
    if ~isempty(requested_checkpoint_dir)
        full_checkpoint_dir = requested_checkpoint_dir;
    elseif abs(relative_uncertainty_level - 5e-3) < 1e-12 && ...
            exist(reuse_full_checkpoint_dir, 'dir')
        full_checkpoint_dir = reuse_full_checkpoint_dir;
    else
        full_checkpoint_dir = fullfile(results_dir, ...
            'full_wcracc_frequency_checkpoints');
    end
    fprintf('Full-WCRACC checkpoint directory: %s\n', full_checkpoint_dir);

    resolved_evaluator = which('evaluate_performance');
    expected_evaluator = fullfile(project_root, 'src', 'evaluations', ...
        'evaluate_performance.m');
    if ~strcmpi(resolved_evaluator, expected_evaluator)
        error('Step2:WrongEvaluatePerformance', ...
            ['MATLAB resolves evaluate_performance to %s, but Step 2 ', ...
             'requires %s.'], resolved_evaluator, expected_evaluator);
    end
    % --- Load Data for the current mode ---
    data_dir = fullfile(project_root, 'data', 'SimulateRIR', current_mode);
    para_file = fullfile(data_dir, 'para.mat');
    % if ~exist(para_file, 'file')
    %     continue; % Skip to the next mode
    % end
    load(para_file); % Loads 'para' setting file

    switch current_mode
        case 'temperature'
            param_vector = para.temperature_vector_celsius;
            design_env_para = ctrl_design_env(1);
        case 'position'
            param_vector = para.perturbation_radii_m;
            design_env_para = ctrl_design_env(2);
        case 'snr'
            param_vector = para.snr_vector;
            design_env_para = ctrl_design_env(3);
    end
    num_levels = length(param_vector); % validation experiment counts

    freq_params = para.freq_params;
    target_freqs = freq_params.target_freqs;
    num_target_freqs = length(target_freqs);

    eval_params = struct();
    if isfield(para, 'eval_params'), eval_params = para.eval_params; end

    % For 'snr' and 'position' modes, tmperature(sound speed) is constant. We pre-compute
    % S_matrix using the base sound speed from the para file.
    % For 'temperature' mode, this will serve as a placeholder and will be
    % overwritten inside the loop.
    geometry_array = load(fullfile(project_root, 'data', ...
        'arrayGeometry', 'array_layout.mat'));
    base_temp = 20; % temperature for case 2&3
    eval_params.S_matrix = precompute_steering_matrix(geometry_array.roomArray.BZ_eval, target_freqs, temp2speed(base_temp));
    % --- Pre-allocate results for the current mode ---
    val_results = struct();
    for k = 1:length(algorithms_to_test)
        for m = 1:length(metrics_to_evaluate)
            val_results.(algorithms_to_test{k}).(metrics_to_evaluate{m}) = nan(num_levels, num_target_freqs);
        end
    end

%% algorithm paras config
    % algorithms_to_test = {'ACC','PM','ACC_PM', 'wcRACC','RPM',...
    % 'RACC_PM_GLS','RACC_PM_Subpro', 'POTDC_RACC'};
    freq_params.virtual_src_idx = 13; % target sound field,PM-related
    % The manuscript experiment uses the plane-wave target throughout.
    desired_data = load(fullfile(project_root, 'data', ...
        'ATF_desired_plane.mat'), 'ATF_desired_plane');
    ATF_desired = desired_data.ATF_desired_plane;


   design_filename = get_data_filename(data_dir, current_mode, design_env_para);
   design_data = load(design_filename);
   ATF_BZ_ctrl = design_data.ATF_BZ.ctrl;
   ATF_DZ_ctrl = design_data.ATF_DZ.ctrl;

   % eta_Z(f)=nu*||H_Z(f)||_F is calculated only from the nominal design
   % ATFs. No off-nominal validation data are used to choose the radius.
   eta_B = local_relative_radii( ...
       ATF_BZ_ctrl, relative_uncertainty_level);
   eta_D = local_relative_radii( ...
       ATF_DZ_ctrl, relative_uncertainty_level);
   freq_params.epsilon = struct('B', eta_B, 'D', eta_D);
   freq_params.gamma = struct('B', eta_B.^2, 'D', eta_D.^2);
   freq_params.scale = 1;
   freq_params.fre_idces = 1:num_target_freqs;
   freq_params.full_racc = struct( ...
       'relative_eta', relative_uncertainty_level, ...
       'diagonal_loading_ratio', diagonal_loading_ratio, ...
       'formulation', 'decomposed', ...
       'solution_method', 'normalized', ...
       'solver', 'mosek', ...
       'mosek_threads', mosek_threads_per_worker, ...
       'parallel_workers', full_wcracc_parallel_workers, ...
       'checkpoint_dir', full_checkpoint_dir, ...
       'resume_from_checkpoints', true, ...
       'retry_unsolved_checkpoints', true, ...
       'normalization_multipliers', [0.1, 1, 5, 10, 0.01], ...
       'fallback_to_bisection', true, ...
       'recovery_method', 'paper', ...
       'randomization_trials', 1000, ...
       'random_seed', 20260820, ...
       'refinement_iterations', 300, ...
       'refinement_seeds', 8);
   freq_params.potdc = struct( ...
       'iterations', 20, ...
       'tolerance', 1e-8, ...
       'solver', 'mosek', ...
       'mosek_threads', mosek_threads_per_worker, ...
       'parallel_workers', full_wcracc_parallel_workers);
   freq_params.racc_pm = struct( ...
       'implementation', 'RACC_PM_Sub', ...
       'parallel_workers', racc_pm_parallel_workers, ...
       'mosek_threads', racc_pm_mosek_threads_per_worker, ...
       'checkpoint_dir', fullfile(results_dir, ...
           'racc_pm_frequency_checkpoints'));

   % ACC_PM:
   kappa = 0.5; % the weight for BZ
   freq_params.kappa = kappa;
   % POTDC_RACC/wcRACC/RPM: bound paras above

   % Use the strict NoCT model as the matched reference for the fixed AC
   % constraint. This avoids importing the historical WCRACC radius.
   noct_reference = NoCT_WCRACC( ...
       ATF_BZ_ctrl, ATF_DZ_ctrl, freq_params);
   filters_w = noct_reference.w;
   alpha_margin_db = 0;
   alpha_AC = zeros(num_target_freqs,1);
   for i = 1:num_target_freqs
       w_f   = filters_w(:, i);
       H_B_f = ATF_BZ_ctrl(:, :, i);
       H_D_f = ATF_DZ_ctrl(:, :, i);
       alpha_AC(i) = 10^((real(calculate_AC( ...
           w_f, H_B_f, H_D_f)) - alpha_margin_db) / 10);
   end
   freq_params.mu = 1; % equal for RPM-ACC & PM
   freq_params.rho = racc_pm_rho;
   freq_params.alpha = alpha_AC;
   para_gamma = freq_params.mu + freq_params.rho * freq_params.alpha;
   freq_params.para_gamma = para_gamma;
   fprintf(['RACC-PM alpha is fixed from nominal control points with ', ...
       'a %.1f dB margin. No evaluation ATF was used.\n'], alpha_margin_db);

   %% --- Get control filters ---
   all_filters = struct();
   algorithm_details = struct();
   algorithm_wall_seconds = struct();
   design_param = design_env_para;
   % load rir/ATF data to design control filter
   design_filename = get_data_filename(data_dir, current_mode, design_param);
   % % test mean value as the data center
   % load('mean_design_data.mat');
   % design_data.ATF_BZ.ctrl = Mean_Data.BZ; design_data.ATF_DZ.ctrl = Mean_Data.DZ;
   %
   design_data.ATF_desired = ATF_desired;
   %
   for k = 1:length(algorithms_to_test)
       algo_name = algorithms_to_test{k};
       fprintf('====== Designing %s (%d/%d) ======\n', ...
           algo_name, k, length(algorithms_to_test));
        design_checkpoint_file = fullfile(design_checkpoint_dir, ...
            sprintf('%02d_%s.mat', k, algo_name));
        restored_design = false;
        if exist(design_checkpoint_file, 'file')
            loaded_checkpoint = load(design_checkpoint_file);
            local_validate_design_checkpoint(loaded_checkpoint, algo_name, ...
                relative_uncertainty_level, eta_B, eta_D, alpha_AC, ...
                alpha_margin_db, freq_params.rho, freq_params.mu);
            all_filters.(algo_name) = loaded_checkpoint.algorithm_filter;
            algorithm_details.(algo_name) = ...
                loaded_checkpoint.algorithm_detail;
            algorithm_wall_seconds.(algo_name) = ...
                loaded_checkpoint.algorithm_wall_time_seconds;
            restored_design = true;
            fprintf('%s restored from algorithm checkpoint.\n', algo_name);
        else
            algorithm_timer = tic;
       % fre_idces = 54; % GLS:7 indices： 8;20;79;104;137;156;158; Sub: 57
       % freq_params.fre_idces = fre_idces;
       if strcmp(algo_name, 'RACC_PM_Subpro')
           algorithm_details.(algo_name) = ...
               run_racc_pm_sub_parallel_checkpointed( ...
               ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, freq_params, ...
               freq_params.racc_pm.checkpoint_dir, ...
               freq_params.racc_pm.parallel_workers);
           all_filters.(algo_name) = algorithm_details.(algo_name).w;
       else
           [all_filters.(algo_name), algorithm_details.(algo_name)] = ...
               design_filters(algo_name, design_data, freq_params);
       end
            algorithm_wall_seconds.(algo_name) = toc(algorithm_timer);
        end
       if strcmp(algo_name, 'Full_WCRACC')
           solved = cellfun(@local_status_is_solved, ...
               algorithm_details.(algo_name).status);
           if any(~solved)
               error('Step2:FullWCRACCUnsolved', ...
                   ['Full_WCRACC was not solved at every frequency. ', ...
                    'Inspect algorithm_details.Full_WCRACC.status ', ...
                    'before deciding whether to use a fallback.']);
           end
       elseif strcmp(algo_name, 'POTDC_RACC')
           solved = cellfun(@local_status_is_solved, ...
               algorithm_details.(algo_name).status);
           if any(~solved)
               error('Step2:POTDCUnsolved', ...
                   ['POTDC_RACC was not solved at every frequency. ', ...
                    'Inspect algorithm_details.POTDC_RACC.status.']);
           end
       elseif strcmp(algo_name, 'RACC_PM_Subpro')
           solved = arrayfun(@(s) local_status_is_solved(char(s)), ...
               algorithm_details.(algo_name).status);
           if any(~solved)
               error('Step2:RACCPMSubUnsolved', ...
                   ['RACC_PM_Subpro was not solved at every frequency. ', ...
                    'Inspect algorithm_details.RACC_PM_Subpro.status.']);
           end
       end
        if ~restored_design
            algorithm_filter = all_filters.(algo_name);
            algorithm_detail = algorithm_details.(algo_name);
            algorithm_wall_time_seconds = ...
                algorithm_wall_seconds.(algo_name);
            checkpoint_rho = freq_params.rho;
            checkpoint_mu = freq_params.mu;
            save(design_checkpoint_file, 'algorithm_filter', ...
                'algorithm_detail', 'algorithm_wall_time_seconds', ...
                'algo_name', 'relative_uncertainty_level', 'eta_B', ...
                'eta_D', 'alpha_AC', 'alpha_margin_db', ...
                'checkpoint_rho', 'checkpoint_mu', '-v7.3');
            fprintf('%s design completed in %.2f minutes.\n', algo_name, ...
                algorithm_wall_seconds.(algo_name) / 60);
        end
    end

 %% Evaluations: simulate the different work condition
 fprintf('======  Evaluating ======\n');
        evaluation_checkpoint_dir = fullfile(results_dir, ...
            'evaluation_checkpoints');
        if ~exist(evaluation_checkpoint_dir, 'dir')
            mkdir(evaluation_checkpoint_dir);
        end
        % using a filter obtained in the nominal environment
        for j = 1:num_levels  % serve for evaluation
            evaluation_checkpoint_file = fullfile( ...
                evaluation_checkpoint_dir, sprintf( ...
                'temperature_%03d.mat', j));
            if exist(evaluation_checkpoint_file, 'file')
                loaded_evaluation = load(evaluation_checkpoint_file);
                if loaded_evaluation.j ~= j || ...
                        ~isequaln(loaded_evaluation.operating_param, ...
                            param_vector(j)) || ...
                        ~isequaln(loaded_evaluation.algorithms_to_test, ...
                            algorithms_to_test) || ...
                        ~isequaln(loaded_evaluation.metrics_to_evaluate, ...
                            metrics_to_evaluate) || ...
                        ~isequaln( ...
                            loaded_evaluation.relative_uncertainty_level, ...
                            relative_uncertainty_level)
                    error('FairEvaluation:EvaluationCheckpointMismatch', ...
                        'Evaluation checkpoint %d has a different setup.', j);
                end
                val_results = loaded_evaluation.val_results;
                fprintf('Evaluation %.1f deg C (%d/%d) restored.\n', ...
                    param_vector(j), j, num_levels);
                continue;
            end
            % log_message(log_fid, sprintf('Processing: Mode [%s], Design [%d/%d], Operating [%d/%d]...', ...
            %     upper(current_mode), i, num_levels, j, num_levels), 'INFO', ECHO_TO_CONSOLE);
            % load cross-evaluation data
            operating_param = param_vector(j); % mode paras values
            operating_filename = get_data_filename(data_dir, current_mode, operating_param);
            operating_data = load(operating_filename); % load the evaluation rir/ATF data
            operating_data.ATF_desired = ATF_desired;

            if strcmp(current_mode, 'temperature')
                    % S_matrix depends on sound speed, which varies with temperature.
                    % Re-calculate it for the current operating condition.
                    current_tmperature = operating_param;
                    eval_params.S_matrix = precompute_steering_matrix...
                        (geometry_array.roomArray.BZ_eval,target_freqs, ...
                        temp2speed(current_tmperature));
            end

            for k = 1:length(algorithms_to_test)
            % for k = 7:7 % for debug
            % algorithms_to_test = {'ACC','PM','ACC_PM', 'wcRACC','RPM',...
            % 'RACC_PM_GLS','RACC_PM_Subpro', 'POTDC_RACC'};
                algo_name = algorithms_to_test{k};
                filters_w = all_filters.(algo_name);
                eval_params.algorithm_name = algo_name;
                performance = evaluate_performance(filters_w, operating_data, freq_params, eval_params);

                for m = 1:length(metrics_to_evaluate)
                    metric_name = metrics_to_evaluate{m};
                    val_results.(algo_name).(metric_name)(j,:) = performance.(metric_name)';
                end
            end
            save(evaluation_checkpoint_file, 'val_results', ...
                'operating_param', 'j', 'algorithms_to_test', ...
                'metrics_to_evaluate', 'relative_uncertainty_level', ...
                '-v7.3');
            fprintf('Evaluated %.1f deg C (%d/%d).\n', ...
                operating_param, j, num_levels);
        end

    % --- Save results for the current mode ---
    filename = "performance_matrices_" + string(run_timestamp) + ".mat";
    results_filename = fullfile(results_dir, filename);
    uncertainty_config = struct( ...
        'relative_level', relative_uncertainty_level, ...
        'eta_B', eta_B, 'eta_D', eta_D, ...
        'gamma_B', eta_B.^2, 'gamma_D', eta_D.^2, ...
        'diagonal_loading_ratio', diagonal_loading_ratio, ...
        'potdc_mapping', ...
        'bright eta_B; dark gamma_D=eta_D^2 (source-paper hybrid model)', ...
        'full_recovery_method', 'paper', ...
        'full_parallel_workers', full_wcracc_parallel_workers, ...
        'racc_pm_parallel_workers', racc_pm_parallel_workers, ...
        'mosek_threads_per_worker', mosek_threads_per_worker, ...
        'full_checkpoint_dir', full_checkpoint_dir, ...
        'racc_pm_checkpoint_dir', freq_params.racc_pm.checkpoint_dir, ...
        'racc_pm_implementation', which('RACC_PM_Sub'), ...
        'algorithm_profile', algorithm_profile, ...
        'alpha', alpha_AC, ...
        'alpha_source', 'nominal NoCT-WCRACC control-point AC', ...
        'alpha_margin_db', alpha_margin_db, ...
        'rho', freq_params.rho, 'mu', freq_params.mu, ...
        'evaluate_performance_path', resolved_evaluator, ...
        'nominal_environment', design_env_para);
    save(results_filename, 'val_results', 'para', 'freq_params', ...
        'algorithms_to_test', 'metrics_to_evaluate', 'all_filters', ...
        'algorithm_details', 'algorithm_wall_seconds', ...
        'uncertainty_config', '-v7.3');
    marker_id = fopen(fullfile(results_dir, 'RUN_COMPLETE.txt'), 'w');
    if marker_id >= 0
        fprintf(marker_id, 'Completed: %s\n', ...
            char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
        fprintf(marker_id, 'Result: %s\n', results_filename);
        fprintf(marker_id, 'RACC-PM: no-suffix RACC_PM_Sub\n');
        fprintf(marker_id, 'relative uncertainty level: %.9g\n', ...
            relative_uncertainty_level);
        fprintf(marker_id, 'rho: %.9g\n', freq_params.rho);
        fprintf(marker_id, 'algorithm profile: %s\n', algorithm_profile);
        fprintf(marker_id, 'alpha margin: %.3f dB\n', alpha_margin_db);
        fclose(marker_id);
    end
end
% --- Final Success Notification ---
elapsed_time_total_hours = toc(total_timer) / 3600;
fprintf('All Validation tasks finished in %.2f hours.\n', elapsed_time_total_hours);
diary off

function radii = local_relative_radii(H, relative_level)
    num_frequencies = size(H, 3);
    radii = zeros(num_frequencies, 1);
    for k = 1:num_frequencies
        radii(k) = relative_level * norm(H(:, :, k), 'fro');
    end
end

function solved = local_status_is_solved(status)
    solved = strcmp(status, 'Solved') || strcmp(status, 'Inaccurate/Solved');
end

function local_validate_design_checkpoint(loaded, algorithm_name, ...
        relative_level, eta_B, eta_D, alpha, alpha_margin_db, rho, mu)
    required = {'algorithm_filter', 'algorithm_detail', ...
        'algorithm_wall_time_seconds', 'algo_name', ...
        'relative_uncertainty_level', 'eta_B', 'eta_D', 'alpha_AC', ...
        'alpha_margin_db', 'checkpoint_rho', 'checkpoint_mu'};
    for k = 1:numel(required)
        if ~isfield(loaded, required{k})
            error('FairEvaluation:IncompleteDesignCheckpoint', ...
                'Design checkpoint for %s lacks %s.', ...
                algorithm_name, required{k});
        end
    end
    if ~strcmp(char(string(loaded.algo_name)), algorithm_name) || ...
            ~isequaln(loaded.relative_uncertainty_level, relative_level) || ...
            ~isequaln(loaded.eta_B, eta_B) || ...
            ~isequaln(loaded.eta_D, eta_D) || ...
            ~isequaln(loaded.alpha_AC, alpha) || ...
            ~isequaln(loaded.alpha_margin_db, alpha_margin_db) || ...
            ~isequaln(loaded.checkpoint_rho, rho) || ...
            ~isequaln(loaded.checkpoint_mu, mu)
        error('FairEvaluation:DesignCheckpointMismatch', ...
            ['The saved %s design does not match the frozen uncertainty ', ...
             'or control-derived alpha configuration.'], algorithm_name);
    end
end

function value = local_environment_scalar(name, default_value)
    text_value = strtrim(getenv(name));
    if isempty(text_value)
        value = default_value;
        return;
    end
    value = str2double(text_value);
    if ~isscalar(value) || ~isfinite(value)
        error('FairEvaluation:InvalidEnvironmentScalar', ...
            '%s must contain one finite numeric scalar.', name);
    end
end
