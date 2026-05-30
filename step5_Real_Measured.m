
% % filter calculations
clear; clc; close all;
addpath(genpath('src/'));
addpath(genpath('data/'));


%% algorithm paras config
algorithms_to_test = {'ACC','PM', 'wcRACC','POTDC_RACC','RPM','RACC_PM_Subpro'}; 
% 'ACC','PM', 'wcRACC','POTDC_RACC','RPM','RACC_PM_Subpro'
% noshow: 'ACC_PM'= 'PM','RACC_PM_GLS',
metrics_to_evaluate = {'AC','NSRE', 'AE'}; % 

total_timer = tic;
% Load Data
data_dir = fullfile('data', 'Cabin_Measurements');
[gamma,epsilon,freq_params] = get_real_measurement_bound(data_dir);
target_freqs = freq_params.target_freqs;
num_levels = freq_params.file_num;
num_target_freqs = length(target_freqs);
virtual_src_idx = 20;
plane_wav_flag = 0;
if plane_wav_flag
    load('./data/ATF_desired_plane.mat'); 
    ATF_desired = ATF_desired_plane;
else
    load('./data/Cabin_Measurements/ATF_2.mat'); % Nominal_ATF.mat
    ATF_desired = squeeze(ATF_BZ(:, virtual_src_idx, :)); 
    % use nominal env
end
% Pre-allocate results 
results = struct();
for k = 1:length(algorithms_to_test)
    for m = 1:length(metrics_to_evaluate)
        results.(algorithms_to_test{k}).(metrics_to_evaluate{m}) = nan(num_levels, num_target_freqs);
    end
end

 % get uncetaity paras for ACC-related algorithms.
 freq_params.epsilon = epsilon;
 freq_params.gamma = gamma;
 freq_params.scale = 1e-3;
 % ACC_PM:
 kappa = 0.5; % the weight for BZ
 freq_params.kappa = kappa;
 % POTDC_RACC/wcRACC/RPM: bound paras above
 % RACC_PM_GLS/RACC_PM_Subpro:
 % alpha for RACC-PM 
 ATF_BZ_ctrl = ATF_BZ;ATF_DZ_ctrl = ATF_DZ;
 wopt = wcACC(ATF_BZ_ctrl, ATF_DZ_ctrl,freq_params);
 if isfield(wopt,'scale')
     scale = wopt.scale;
     filters_w = wopt.w;
 else
     scale = nan;
     filters_w = wopt;
 end
 margin = 5; % Leaving a 3 dB margin to the inaccuracies upper bound from wcRACC.
 alpha_AC = zeros(num_target_freqs,1);
 for i = 1:num_target_freqs
 % Extract frequency-dependent data for this iteration.
     w_f   = filters_w(:, i);
     H_B_f = ATF_BZ(:, :, i); 
     H_D_f = ATF_DZ(:, :, i);
     alpha_AC(i) = 10^((real(calculate_AC(w_f, H_B_f, H_D_f))-margin)/10);
 end
 freq_params.mu = 0.5; % equal for RPM-ACC & PM
 freq_params.rho = 10; % trade-off weights:1000-->ACC, 1e-3-->PM
 freq_params.alpha = alpha_AC; % refer to wcRACC-3dB
 para_gamma = freq_params.mu + freq_params.rho * freq_params.alpha;
 freq_params.para_gamma = para_gamma;
 

   %% --- Get control filters ---
   all_filters = struct();
   % data to design control filter 
   design_data.ATF_BZ.ctrl = ATF_BZ_ctrl;
   design_data.ATF_DZ.ctrl = ATF_DZ_ctrl;
   design_data.ATF_desired = ATF_desired;

   % for k = 1:length(algorithms_to_test) % 7:7 debug
   for k = 1:length(algorithms_to_test)
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
% load('real_measured_filter.mat');
evaluating_data.ATF_BZ = struct(); evaluating_data.ATF_DZ = struct();
evaluating_data.ATF_desired = ATF_desired;
evaluating_data.ATF_BZ.ctrl = ATF_BZ_ctrl;
fprintf('======  Evaluating ======\n');
% using a filter obtained in the nominal environment
for j = 1:num_levels  % serve for evaluation
    operating_filename = sprintf('ATF_%d.mat', j);
    tmp_evaluating_data = load(operating_filename); % load the ATF data
    evaluating_data.ATF_BZ.eval = tmp_evaluating_data.ATF_BZ;
    evaluating_data.ATF_DZ.eval = tmp_evaluating_data.ATF_DZ;            
    % for k = 1:length(algorithms_to_test)
    for k = 1:length(algorithms_to_test)
        algo_name = algorithms_to_test{k};
        filters_w = all_filters.(algo_name);
        eval_params.algorithm_name = algo_name;
        eval_params.virtual_source_idx = virtual_src_idx;
        performance = evaluate_performance(filters_w, evaluating_data, freq_params, eval_params);
        for m = 1:length(metrics_to_evaluate)
            metric_name = metrics_to_evaluate{m};
            val_results.(algo_name).(metric_name)(j,:) = performance.(metric_name)';
        end
    end
end
% --- Save results for the current mode ---
results_dir = fullfile('real_measured_results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
% timestamp_suffix 
timestamp_str = datetime('now', 'format','yyyyMMdd_HHmm');    
% joint filename
filename = "performance_matrices_" + string(timestamp_str) + ".mat";
results_filename = fullfile(results_dir, filename);
save(results_filename, 'val_results', 'algorithms_to_test', 'metrics_to_evaluate','all_filters');
% --- Final Success Notification ---
elapsed_time_total_hours = toc(total_timer) / 3600;
fprintf('All Validation tasks finished in %.2f hours.\n', elapsed_time_total_hours);

%% Plot Settings
% 将算法分为两组
metrics_to_evaluate = {'AC','NSRE'};
group1_algs = {'ACC', 'wcRACC', 'POTDC_RACC', 'RACC_PM_Subpro','PM', 'RPM'};
group2_algs = {'PM', 'RPM'};
all_groups = {group1_algs, group2_algs};
group_names = {'Proposed RACC_PM_Subpro & Others', 'Proposed RPM & Others'};
design_env_idx = 2;
% x-axis
freqs = freq_params.target_freqs;
freq_params.target_freq_end = 5000;

base_colors = [0.00,0.44,0.72;
    0.98,0.69,0.10;
    0.00,0.70,0.69;
    0.83,0.05,0.55;
    0.76, 0.87, 0.69;
    0.76, 0.87, 0.69];
% Marker
markers = {'o', 's', 'd','none','none','>'}; line_style ={'-.',':','-','-','--','-'};

for j = 1:length(metrics_to_evaluate)
    metric_name = metrics_to_evaluate{j};    
    % Group 1  Group 2 figures
    for g = 1:2
        current_algs = all_groups{g};
        
        figure('Name', sprintf('%s Performance - %s', metric_name, group_names{g}), ...
               'Position',[390 + g*50, 560 + g*50, 1400, 480]);
        hold on; grid on;
        % =========================================================
        % mean results
        % =========================================================
        for i = 1:length(current_algs)
            alg_name = current_algs{i};        
            y_data_all = val_results.(alg_name).(metric_name);
            num_total_envs = size(y_data_all, 1);
            val_indices = setdiff(1:num_total_envs, design_env_idx);        
            y_data_val = y_data_all(val_indices, :);
            y_mean = mean(y_data_all, 1);            
            % highlightened mean results lines
            idx_markers = [1 3 6 8 13 20 25 30 50]; 
            semilogx(freqs, y_mean, 'LineStyle',line_style{i}, 'Marker',markers{i}, ...
                'Color', base_colors(i, :), ...
                'LineWidth', 2, ...                
                'MarkerSize', 8, ...               
                'MarkerIndices', idx_markers, ...           
                'DisplayName', strrep(alg_name, '_', '\_'));
            legend('Location', 'best', 'FontSize', 14, 'NumColumns', 4,'FontName','Times New Roman');    
        end  
        % =========================================================
        % error bar
        % =========================================================
        target_sparse_freqs =[500, 1000, 2000, 4000];         
        % for clarity
        shift_factors =[0.92 0.95, 0.98, 1.01, 1.04, 1.07]; 
        
        for i = 1:length(current_algs)
            alg_name = current_algs{i};        
            y_data_all = val_results.(alg_name).(metric_name);
            val_indices = setdiff(1:size(y_data_all, 1), design_env_idx);
            y_data_val = y_data_all(val_indices, :);            
            for f_idx = 1:length(target_sparse_freqs)
           
                [~, actual_idx] = min(abs(freqs - target_sparse_freqs(f_idx)));
                actual_f = freqs(actual_idx);
                
                data_at_f = y_data_val(:, actual_idx);
                
                m_val = mean(data_at_f);
                u_val = prctile(data_at_f, 90);
                l_val = prctile(data_at_f, 10);
                
                % errorbar
                err_neg = m_val - l_val;
                err_pos = u_val - m_val;
                
                shifted_f = actual_f * shift_factors(i);
                
                % errorbar legend
                % if f_idx ==1 
                %     errorbar(shifted_f, m_val, err_neg, err_pos, ...
                %         'LineStyle', 'none', 'Color', base_colors(i, :), ...
                %         'LineWidth', 1.5, 'CapSize', 7, ... 
                %         'HandleVisibility', 'off','Marker',markers{i});
                % else
                    errorbar(shifted_f, m_val, err_neg, err_pos, ...
                        'LineStyle', 'none', 'Color', base_colors(i, :), ...
                        'LineWidth', 1.5, 'CapSize', 7, ... 
                        'HandleVisibility', 'off',...
                        'Marker',markers{i});
                % end
            end
        end
      h_dummy_err = errorbar(10, -100, 5, 5,'Color',[0.3 0.3 0.3], ...
    'LineWidth', 1.5, 'CapSize', 7, 'Marker', 'none','LineStyle','none');
        % =========================================================
        % decorations
        % =========================================================
        xlabel('Frequency: [Hz]', 'FontSize', 15.4,'FontName','Times New Roman');
        ylabel([metric_name ': [dB]'], 'FontSize', 15.4,'FontName','Times New Roman');    
        xlim([freq_params.target_freq_start, freq_params.target_freq_end]);
        
        set(gca, 'XScale','log', 'Box', 'on');
        hold off;
    end
end