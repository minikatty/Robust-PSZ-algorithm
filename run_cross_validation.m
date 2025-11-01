%% MAIN SCRIPT - Automated Cross-Validation for Server
% This script loads a pre-generated RIR/ATF database, runs a full N x N 
% cross-validation for multiple algorithms, and saves the resulting 
% performance matrices. Includes logging and email notifications.
%
% Author: ZhouLei
% Date: 2025/10/27 (Refactored for server execution)

clear; clc; close all;

%% ========================================================================
%  1. CONFIGURATION & SETUP
%  ========================================================================
addpath(genpath('src/'));
rng(2025);

% --- Logging & Notification ---
log_dir = 'logs';
if ~exist(log_dir, 'dir'), mkdir(log_dir); end
run_timestamp = datetime("now","Format","yyyyMMdd_HHmm");
log_filename = fullfile(log_dir, sprintf('%s_cross_validation_run.log', char(run_timestamp)));
log_fid = fopen(log_filename, 'a');
if log_fid == -1, error('main_2:LogOpenFailed', '无法打开日志文件: %s', log_filename); end
cleanupObj = onCleanup(@() safe_close(log_fid));
ECHO_TO_CONSOLE = true; % Set true for local debug, false for server nohup

recipient_email = 'zhouleicqupt2016@outlook.com'; % <-- Set your email
enable_notifications = true;

% --- Experiment Configuration ---
experiment_modes = {'snr'}; %, 'temperature', 'position'
algorithms_to_test = {'ACC'};% , 'PM','ACC-PM', 'WCRACC','RPM', 'RACC-PM', 'POTDC_RACC'
metrics_to_evaluate = {'AC'}; % , 'NSRE', 'AE','Planarity'

%% ========================================================================
%  2. MAIN EXECUTION BLOCK
%  ========================================================================
try
    log_message(log_fid, '--- Cross-Validation Script Started ---', 'INFO', ECHO_TO_CONSOLE);

    total_timer = tic;

    % --- Loop through each experiment mode ---
    for mode_idx = 1:length(experiment_modes)
        current_mode = experiment_modes{mode_idx};
        mode_timer = tic;
        
        log_message(log_fid, sprintf('====== Starting Experiment Mode: %s ======', upper(current_mode)), 'INFO', ECHO_TO_CONSOLE);

        % --- Load Data for the current mode ---
        data_dir = fullfile('data', 'SimulateRIR', current_mode);
        para_file = fullfile(data_dir, 'para.mat');
        if ~exist(para_file, 'file')
            log_message(log_fid, sprintf('Metadata file "para.mat" not found in %s. Skipping this mode.', data_dir), 'WARN', ECHO_TO_CONSOLE);
            continue; % Skip to the next mode
        end
        load(para_file); % Loads 'para'
        
        switch current_mode
            case 'snr', param_vector = para.snr_vector;
            case 'temperature', param_vector = para.temperature_vector_celsius;
            case 'position', param_vector = para.perturbation_radii_m;
        end
        num_levels = length(param_vector);
        freq_params = para.freq_params;

        % --- Pre-allocate results for the current mode ---
        results = struct();
        for k = 1:length(algorithms_to_test)
            for m = 1:length(metrics_to_evaluate)
                results.(algorithms_to_test{k}).(metrics_to_evaluate{m}) = nan(num_levels, num_levels);
            end
        end
        log_message(log_fid, sprintf('Initialized %dx%d results matrices for %d algorithms.', ...
            num_levels, num_levels, length(algorithms_to_test)), 'INFO', ECHO_TO_CONSOLE);

        % --- Main Cross-Validation Loop ---
        for i = 1:num_levels
            design_param = param_vector(i);
            design_filename = get_data_filename(data_dir, current_mode, design_param);
            design_data = load(design_filename);
            
            all_filters = struct();
            for k = 1:length(algorithms_to_test)
                algo_name = algorithms_to_test{k};
                all_filters.(algo_name) = design_filters(algo_name, design_data, freq_params);
            end
            
            % cross-evaluations
            for j = 1:num_levels
                log_message(log_fid, sprintf('Processing: Mode [%s], Design [%d/%d], Operating [%d/%d]...', ...
                    upper(current_mode), i, num_levels, j, num_levels), 'INFO', ECHO_TO_CONSOLE);
                 % load cross-evaluation data
                if i == j, operating_data = design_data;
                else
                    operating_param = param_vector(j); % snr sets
                    operating_filename = get_data_filename(data_dir, current_mode, operating_param);
                    operating_data = load(operating_filename);
                end
                
                for k = 1:length(algorithms_to_test)
                    algo_name = algorithms_to_test{k};
                    filters_w = all_filters.(algo_name);
                    performance = evaluate_performance(filters_w, operating_data, freq_params);
                    
                    for m = 1:length(metrics_to_evaluate)
                        metric_name = metrics_to_evaluate{m};
                        results.(algo_name).(metric_name)(i, j) = performance.(metric_name);
                    end
                end
            end
        end
        
        % --- Save results for the current mode ---
        results_dir = fullfile('results', [current_mode, '_AnalysisData']);
        if ~exist(results_dir, 'dir'), mkdir(results_dir); end
        results_filename = fullfile(results_dir, 'performance_matrices.mat');
        save(results_filename, 'results', 'para', 'algorithms_to_test', 'metrics_to_evaluate');
        
        elapsed_time_mode = toc(mode_timer);
        log_message(log_fid, sprintf('SUCCESS: Mode "%s" completed in %.2f minutes. Results saved to %s', upper(current_mode), elapsed_time_mode/60, results_filename), 'INFO', ECHO_TO_CONSOLE);
    end

    % --- Final Success Notification ---
    elapsed_time_total_hours = toc(total_timer) / 3600;
    final_message = sprintf('All cross-validation tasks finished in %.2f hours.', elapsed_time_total_hours);
    log_message(log_fid, final_message, 'INFO', ECHO_TO_CONSOLE);
    
    if enable_notifications
        send_graphmail(recipient_email, '[MATLAB Job Finished] Cross-Validation Success', final_message, {log_filename});
    end
    
    if ECHO_TO_CONSOLE, fprintf('\n--- All Tasks Finished Successfully ---\n'); end

catch ME_main
    % --- Global Error Handling & Notification ---
    final_error_msg = sprintf('CRITICAL FAILURE in cross-validation script.\n\nError: %s', ME_main.message);
    log_message(log_fid, final_error_msg, 'FATAL', ECHO_TO_CONSOLE);
    log_message(log_fid, getReport(ME_main, 'extended', 'hyperlinks', 'off'), 'DEBUG', ECHO_TO_CONSOLE);
    
    if enable_notifications
        send_graphmail(recipient_email, '[MATLAB Job FAILED] Cross-Validation Error', final_error_msg, {log_filename});
    end
    
    if ECHO_TO_CONSOLE, fprintf('\n--- Script Terminated Due to a Critical Error ---\n'); end
    rethrow(ME_main);
end

%% LOCAL FUNCTIONS
function safe_close(fid)
    if ~(isscalar(fid) && isnumeric(fid) && fid == -1), try, fclose(fid); catch, end, end
end