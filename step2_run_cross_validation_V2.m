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
addpath(genpath('data/'));
% --- Experiment Configuration ---
experiment_modes = {'temperature'}; %, , 'position', 'snr', 'temperature'
% algorithms_to_test = {'ACC'};
algorithms_to_test = {'RPM'}; % 'RACC_PM_GLS',
%'ACC','PM','ACC_PM', 'wcRACC','RPM',...
    % 'RACC_PM_GLS','RACC_PM_Subpro', 'POTDC_RACC'
%{'ACC','PM','ACC-PM', 'wcRACC','RPM', 'RACC-PM', 'POTDC-RACC'};% , 
metrics_to_evaluate = {'AC','NSRE', 'AE','Planarity'}; % , 
ctrl_design_env = [22.5, 0.1, 20]; % the enviroment for each case

%% ========================================================================
%  2. MAIN EXECUTION BLOCK
%  ========================================================================
total_timer = tic;
% --- Loop through each experiment mode ---
for mode_idx = 1:length(experiment_modes)
    current_mode = experiment_modes{mode_idx};
    mode_timer = tic;
    
    fprintf('====== Starting Experiment Mode: %s ======\n', upper(current_mode));
    % --- Load Data for the current mode ---
    data_dir = fullfile('data', 'SimulateRIR', current_mode);
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
    geometry_array = load('data/arrayGeometry/array_layout.mat');
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
    % plane wave or monopole wave
    plane_wav_flag = 1;
    if plane_wav_flag
        load('./data/ATF_desired_plane.mat'); 
        ATF_desired = ATF_desired_plane;
    else
        ATF_desired = squeeze(ATF_BZ_ctrl(:, virtual_src_idx, :)); 
        % use nominal env
    end
   

   design_filename = get_data_filename(data_dir, current_mode, design_env_para);
   % get uncetaity paras for ACC-related algorithms. % potential risk for
   % data leakage
   [gamma, epsilon] = get_bound_paras(current_mode, design_filename);
   % freq_params.epsilon = epsilon;
   % freq_params.gamma = gamma;

   % ACC_PM:
   kappa = 0.5; % the weight for BZ
   freq_params.kappa = kappa;
   % POTDC_RACC/wcRACC/RPM: bound paras above

   % RACC_PM_GLS/RACC_PM_Subpro:
   % alpha for RACC-PM 
   load(design_filename);
   ATF_BZ_ctrl = ATF_BZ.ctrl;ATF_DZ_ctrl = ATF_DZ.ctrl;
   ATF_BZ_eval = ATF_BZ.eval;ATF_DZ_eval = ATF_DZ.eval; %eval,ctrl
   freq_params.scale = 1e-2;
   wopt = wcACC(ATF_BZ_ctrl, ATF_DZ_ctrl,freq_params);
   if isfield(wopt,'scale')
       scale = wopt.scale;
       filters_w = wopt.w;
   else
       scale = nan;
       filters_w = wopt;
   end
   margin = 2; % Leaving a 3 dB margin to the inaccuracies upper bound from wcRACC.
   alpha_AC = zeros(num_target_freqs,1);
   for i = 1:num_target_freqs
   % Extract frequency-dependent data for this iteration.
       w_f   = filters_w(:, i);
       H_B_f = ATF_BZ_eval(:, :, i); 
       H_D_f = ATF_DZ_eval(:, :, i);
       alpha_AC(i) = 10^((real(calculate_AC(w_f, H_B_f, H_D_f))-margin)/10);
   end
   freq_params.mu = 1; % equal for RPM-ACC & PM
   freq_params.rho = 1e-4; % trade-off weights:1000-->ACC, 1e-3-->PM
   freq_params.alpha = alpha_AC; % refer to wcRACC-3dB
   para_gamma = freq_params.mu + freq_params.rho * freq_params.alpha;
   freq_params.para_gamma = para_gamma;
   
   %% --- Get control filters ---
   all_filters = struct();
   design_param = design_env_para;
   % load rir/ATF data to design control filter
   design_filename = get_data_filename(data_dir, current_mode, design_param);
   design_data = load(design_filename);
   % % test mean value as the data center
   % load('mean_design_data.mat');
   % design_data.ATF_BZ.ctrl = Mean_Data.BZ; design_data.ATF_DZ.ctrl = Mean_Data.DZ;
   %
   design_data.ATF_desired = ATF_desired;
   %
   for k = 1:length(algorithms_to_test) 
   % for k = 7:7 % debug
       algo_name = algorithms_to_test{k};
       % fre_idces = 54; % GLS:7 indices： 8;20;79;104;137;156;158; Sub: 57
       % freq_params.fre_idces = fre_idces;
       all_filters.(algo_name) = design_filters(algo_name, design_data, freq_params);
       if isfield(wopt,'scale')
           scale = wopt.scale;
           filters_w = wopt.w;
       else
           scale = nan;
           filters_w = wopt;
       end
       % get control filter in a nominal environment
   end
        
%% Evaluations: simulate the different work condition
fprintf('======  Evaluating ======\n');
        % using a filter obtained in the nominal environment
        for j = 1:num_levels  % serve for evaluation
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
        end
        
    % --- Save results for the current mode ---
    results_dir = fullfile('results', [current_mode, '_AnalysisData']);
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    % timestamp_suffix 
    timestamp_str = datetime('now', 'format','yyyyMMdd_HHmm');    
    % joint filename
    filename = "performance_matrices_" + string(timestamp_str) + ".mat";
    results_filename = fullfile(results_dir, filename);
    save(results_filename, 'val_results', 'para', 'algorithms_to_test', 'metrics_to_evaluate','all_filters');
end
% --- Final Success Notification ---
elapsed_time_total_hours = toc(total_timer) / 3600;
fprintf('All Validation tasks finished in %.2f hours.\n', elapsed_time_total_hours);