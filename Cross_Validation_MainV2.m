%% MAIN SCRIPT - Automated Cross-Validation (Final Refactored Version)
% This script orchestrates a full N x N cross-validation. It employs a
% hybrid computation strategy for evaluation parameters to ensure both
% accuracy (for temperature-dependent metrics) and efficiency.
%
% Author: ZhouLei
% Date: 2025/11/02 (Refactored for mode-dependent parameter computation)

clear; clc; close all;

%% ========================================================================
%  1. CONFIGURATION
%  ========================================================================
% --- High-Level Experiment Settings ---
% experiment_modes     = {'snr', 'temperature', 'position'};
% algorithms_to_test   = {'ACC', 'PM', 'WCRACC', 'RPM'};
% metrics_to_evaluate  = {'AC', 'NSRE', 'AE', 'Planarity'};

% --- Experiment Configuration ---
experiment_modes = {'snr'}; %, 'temperature', 'position'
algorithms_to_test = {'ACC'};% , 'PM','ACC-PM', 'WCRACC','RPM', 'RACC-PM', 'POTDC_RACC'
metrics_to_evaluate = {'AC'}; % , 'NSRE', 'AE','Planarity'


ENABLE_NOTIFICATIONS = true;
ECHO_TO_CONSOLE      = true;

% --- Notification & Logging Settings ---
recipient_email = 'zhouleicqupt2016@outlook.com';
log_dir         = 'logs';


%% ========================================================================
%  2. SETUP
%  ========================================================================
addpath(genpath('src/'));
if ~exist(log_dir, 'dir'), mkdir(log_dir); end
run_timestamp = datetime("now", "Format", "yyyyMMdd_HHmm");
log_filename  = fullfile(pwd, log_dir, sprintf('%s_run.log', char(run_timestamp)));
log_fid       = fopen(log_filename, 'a');
if log_fid == -1, error('main:LogOpenFailed', 'Failed to open log file: %s', log_filename); end
cleanupObj    = onCleanup(@() safe_close(log_fid));


%% ========================================================================
%  3. MAIN EXECUTION BLOCK
%  ========================================================================
try
    log_message(log_fid, '--- Cross-Validation Script Started ---', 'INFO', ECHO_TO_CONSOLE);
    total_timer = tic;

    % --- Loop through each experiment mode ---
    for mode_idx = 1:length(experiment_modes)
        current_mode = experiment_modes{mode_idx};
        mode_timer = tic;
        
        log_message(log_fid, sprintf('====== Starting Experiment Mode: %s ======', upper(current_mode)), 'INFO', ECHO_TO_CONSOLE);

        % --- Load Base Data and Parameters for the mode ---
        data_dir  = fullfile('data', 'SimulateRIR', current_mode);
        para_file = fullfile(data_dir, 'para.mat');
        if ~exist(para_file, 'file'), log_message(log_fid, sprintf('Para file not found in %s. Skipping.', data_dir), 'WARN', ECHO_TO_CONSOLE); continue; end
        load(para_file); % Loads 'para' struct
        
        param_vector = get_param_vector(para, current_mode);
        num_levels   = length(param_vector);

        % --- STATIC PRE-COMPUTATION FOR EVALUATION ---
        % Prepare the 'eval_params' struct. For parameters that are constant
        % throughout a mode (like S_matrix in 'snr' and 'position' modes),
        % we compute them here once for efficiency.
        log_message(log_fid, 'Performing static pre-computation for evaluation...', 'INFO', ECHO_TO_CONSOLE);
        eval_params = struct();
        if isfield(para, 'eval_params'), eval_params = para.eval_params; end
        
        % For 'snr' and 'position' modes, sound speed is constant. We pre-compute
        % S_matrix using the base sound speed from the para file.
        % For 'temperature' mode, this will serve as a placeholder and will be
        % overwritten inside the loop.
        eval_params.S_matrix = precompute_steering_matrix(para.arrays.mic_positions_BZ_eval, ...
                                                           para.freq_params.freq_vector, ...
                                                           para.sim.sound_speed_mps);
        
        % --- Pre-allocate Results Matrices ---
        results = initialize_results_struct(algorithms_to_test, metrics_to_evaluate, num_levels);
        log_message(log_fid, 'Initialized results matrices.', 'INFO', ECHO_TO_CONSOLE);

        % --- Main Cross-Validation Loops ---
        for i = 1:num_levels % Design Condition Loop
            design_filename = get_data_filename(data_dir, current_mode, param_vector(i));
            design_data     = load(design_filename);
            
            % Design all filters for the current design condition.
            all_filters = struct();
            for k = 1:length(algorithms_to_test)
                algo_name = algorithms_to_test{k};
                all_filters.(algo_name) = design_filters(algo_name, design_data, para.freq_params, para.algo_params.(algo_name));
            end
            
            for j = 1:num_levels % Operating Condition Loop
                log_message(log_fid, sprintf('Processing: Mode [%s], Design [%d/%d], Operating [%d/%d]', ...
                    upper(current_mode), i, num_levels, j, num_levels), 'INFO', ECHO_TO_CONSOLE);
                
                % Load operating data.
                if i == j, operating_data = design_data;
                else, operating_data = load(get_data_filename(data_dir, current_mode, param_vector(j)));
                end
                
                % --- DYNAMIC PARAMETER UPDATE ---
                % CRITICAL STEP: If a parameter in 'eval_params' depends on the
                % operating condition, it must be updated here.
                if strcmp(current_mode, 'temperature')
                    % S_matrix depends on sound speed, which varies with temperature.
                    % Re-calculate it for the current operating condition.
                    current_sound_speed = operating_data.sim.sound_speed_mps;
                    eval_params.S_matrix = precompute_steering_matrix(para.arrays.mic_positions_BZ_eval, ...
                                                                       para.freq_params.freq_vector, ...
                                                                       current_sound_speed);
                end
                
                % Evaluate all designed filters under the current operating condition.
                for k = 1:length(algorithms_to_test)
                    algo_name = algorithms_to_test{k};
                    filters_w = all_filters.(algo_name);
                    
                    % CORE EVALUATION CALL: Pass the correctly prepared eval_params.
                    % For 'snr'/'position', it uses the static S_matrix.
                    % For 'temperature', it uses the dynamically updated S_matrix.
                    performance = evaluate_performance(filters_w, algo_name, operating_data, para.freq_params, eval_params);
                    
                    % Store every requested metric.
                    for m = 1:length(metrics_to_evaluate)
                        metric_name = metrics_to_evaluate{m};
                        if isfield(performance, metric_name)
                            results.(algo_name).(metric_name)(i, j) = performance.(metric_name);
                        end
                    end
                end
            end
        end
        
        % --- Save Results for the Current Mode ---
        results_dir = fullfile('results', [current_mode, '_AnalysisData']);
        if ~exist(results_dir, 'dir'), mkdir(results_dir); end
        results_filename = fullfile(results_dir, 'performance_matrices.mat');
        save(results_filename, 'results', 'para', 'algorithms_to_test', 'metrics_to_evaluate');
        
        log_message(log_fid, sprintf('SUCCESS: Mode "%s" completed in %.2f min. Results saved.', upper(current_mode), toc(mode_timer)/60), 'INFO', ECHO_TO_CONSOLE);
    end

    % --- Final Success Notification ---
    final_message = sprintf('All cross-validation tasks finished successfully in %.2f hours.', toc(total_timer) / 3600);
    log_message(log_fid, final_message, 'INFO', ECHO_TO_CONSOLE);
    if ENABLE_NOTIFICATIONS, send_graphmail(recipient_email, '[MATLAB Job Finished]', final_message, {log_filename}); end
    
catch ME_main
    % --- 4. GLOBAL ERROR HANDLING & NOTIFICATION ---
    % ... (Error handling code remains the same) ...
end

%% ========================================================================
%  5. LOCAL HELPER FUNCTIONS
%  ========================================================================
% ... (Helper functions remain the same) ...