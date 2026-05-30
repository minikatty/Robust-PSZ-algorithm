clear; clc; close all;
% experiment_modes = {'snr'}; %, , 'position', 'snr', 'temperature'
ctrl_design_env_idx = [26, 26, 26]; % corresponding to the env [22.5, 0.1, 20]; temperature/positon/snr
eval_env_idx = [18,1,1];% corresponding to the env [25, ~, ~]; temperature/positon/snr
eval_env_ATF_data = "D:\CodeManage\Robust-PSZ-algorithm_V2\data\SimulateRIR\temperature\Data_T-20.00.mat";
current_mode = 'temperature';
ATFs_monitor_file ="D:\CodeManage\Robust-PSZ-algorithm_V2\data\MonitorGrid\Monitor_ATF.mat";
f_target = 1000;
switch current_mode
    case 'temperature'
        % param_vector = para.temperature_vector_celsius;
        design_env_idx = ctrl_design_env_idx(1);
        eval_work_env_idx = eval_env_idx(1);
    case 'position'
        % param_vector = para.perturbation_radii_m;
        design_env_idx = ctrl_design_env_idx(2);
        eval_work_env_idx = eval_env_idx(2);
    case 'snr' 
        % param_vector = para.snr_vector;
        design_env_idx = ctrl_design_env_idx(3);
        eval_work_env_idx = eval_env_idx(3);
end

addpath(genpath('results/temperature_AnalysisData'));
load('performance_matrices_20260430_2119.mat');
performance = val_results;
algorithms_to_test = {'PM','RPM'};
% algorithms_to_test = { 'ACC','wcRACC','POTDC_RACC',...
%     'RACC_PM_Subpro','PM', 'RPM'};
% algorithms_to_test = {'RACC_PM_Subpro'};
% 'ACC','wcRACC','POTDC_RACC',...
    % 'RACC_PM_Subpro','PM', 'RPM'
metrics_to_evaluate = {'Planarity'}; % ,'NSRE', 'AE','Planarity','AC','AE','Planarity'
% some common paras
load("D:\CodeManage\Robust-PSZ-algorithm_V2\data\arrayGeometry\array_layout.mat");
s = roomArray.s;
freq_params = para.freq_params;
rp = 0.8; c = 20;
f_max_val = temp2speed(c)*size(s,1)/(4*pi*rp);

base_colors = [0.00,0.44,0.72;
    0.98,0.69,0.10;
    0.00,0.70,0.69;
    0.83,0.05,0.55;
    0.76, 0.87, 0.69;
    0.76, 0.87, 0.69];
% algorithms_to_test = {'ACC', 'wcRACC','POTDC_RACC','RPM','RACC_PM_Subpro'}; 
% real_measured algorithms_to_test = {'ACC','wcRACC','POTDC_RACC', 'RACC_PM_Subpro','PM', 'RPM'};
% Marker
markers = {'o', 's', 'd','none','none','>'}; 
line_style ={'-.',':','-','-','--','-'};
idx_markers = 9:10:159; 
idx_markers =[1 idx_markers];
freq_params.target_freq_end = 8000;

for j = 1:length(metrics_to_evaluate)
    metric_name = metrics_to_evaluate{j}; % 
    figure('Name', [metric_name, ' Performance'], 'Position',[100, 100, 800, 500]);
    hold on; grid on;
    if strcmp(metric_name,'NSRE')
        base_colors2 = [0.83,0.05,0.55;
        0.76, 0.87, 0.69;
        0.76, 0.87, 0.69];
        algorithms_to_test_PMtype = {'PM', 'RPM'}; % 'PM', 'RPM'
        for i = 1:length(algorithms_to_test_PMtype)
            alg_name = algorithms_to_test_PMtype{i};        
            y_data = performance.(alg_name).(metric_name)(design_env_idx,:);        % (ctrl_design_env_idx)
            semilogx(freq_params.target_freqs, y_data, ...
                'Color', base_colors2(i, :), ...
                'Marker',markers{i}, ...
                'LineStyle',line_style{i}, ...
                'LineWidth', 1.5, 'MarkerSize', 6, ...
                'MarkerIndices', idx_markers, ... 
                'DisplayName', strrep(alg_name, '_', '-')); % strrep 防止下划线变成下标
            hold on;
             xl = xline(f_max_val, '--r', '$f_{max}\approx 1637$',...
            'LineWidth', 1.5, 'FontName','Times New Roman',... 
            'FontSize', 10, ...
            'LabelVerticalAlignment', 'bottom', ...
            'LabelHorizontalAlignment', 'center',...
            'LabelOrientation', 'horizontal'); % 
            xl.Annotation.LegendInformation.IconDisplayStyle = 'off';
            xl.Interpreter = 'latex';
        end
    else
        for i = 1:length(algorithms_to_test)
            alg_name = algorithms_to_test{i};        
            y_data = performance.(alg_name).(metric_name)(design_env_idx,:);        % (ctrl_design_env_idx)
            semilogx(freq_params.target_freqs, y_data, ...
                'Color', base_colors(i, :), ...
                'Marker',markers{i}, ...
                'LineStyle',line_style{i}, ...
                'LineWidth', 1.5, 'MarkerSize', 6, ...
                'MarkerIndices', idx_markers, ... 
                'DisplayName', strrep(alg_name, '_', '-')); % strrep 防止下划线变成下标
            hold on;
             xl = xline(f_max_val, '--r', '$f_{max}\approx 1637$',...
            'LineWidth', 1.5, 'FontName','Times New Roman',... 
            'FontSize', 10, ...
            'LabelVerticalAlignment', 'bottom', ...
            'LabelHorizontalAlignment', 'center',...
            'LabelOrientation', 'horizontal'); % 
            xl.Annotation.LegendInformation.IconDisplayStyle = 'off';
            xl.Interpreter = 'latex';
        end
    end
    title(sprintf('%s Performance Comparison', metric_name), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Target Frequency (Hz)', 'FontSize', 12);
    ylabel(metric_name, 'FontSize', 12);    
    xlim([freq_params.target_freq_start  freq_params.target_freq_end]);

    legend('Location', 'best', 'FontSize', 10, 'NumColumns', 2);
    hold off;
    figure('Name', [metric_name, ' Performance'], 'Position',[100, 100, 800, 500]);
    hold on; grid on;

    for i = 1:length(algorithms_to_test)
        alg_name = algorithms_to_test{i};        
        y_data = performance.(alg_name).(metric_name)(eval_work_env_idx,:);        % (eval_work_env_idx)        
        semilogx(freq_params.target_freqs, y_data, 'Marker', markers{i}, ...
            'LineWidth', 1.5, 'MarkerSize', 6, ...
            'MarkerIndices', idx_markers, ... 
            'DisplayName', strrep(alg_name, '_', '-')); % strrep 防止下划线变成下标
    end

    title(sprintf('%s Performance Comparison (Validation)', metric_name), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Target Frequency (Hz)', 'FontSize', 12);
    ylabel(metric_name, 'FontSize', 12);
    xlim([freq_params.target_freq_start freq_params.target_freq_end]);
    legend('Location', 'best', 'FontSize', 10, 'NumColumns', 2);
end
%% plot sound field distribution
% for i = 7:7%1:length(algorithms_to_test)
%     % length(algorithms_to_test)
%     alg_name = algorithms_to_test{i};
%     filters_w = all_filters.(alg_name);
%     plot_soundfield_pressure(ATFs_monitor_file,eval_env_ATF_data, filters_w,f_target, ...
%     'SpeakerDataFile',"data/arrayGeometry/array_layout.mat");
%     hold off;
% end

%% plot metrics

% ylim([freq_params.target_freq_start freq_params.target_freq_end]);
% 
% all_axes = findall(gcf, 'type', 'axes'); 
% set(all_axes, ...
%     'XScale', 'log', ...              
%     'XGrid', 'on', ...                
%     'XMinorGrid', 'on', ...           
%     'YGrid', 'on', ...
%     'FontName', 'Times New Roman',...
%     'FontSize', 10, ...               
%     'LineWidth', 1, ...               
%     'XTick', [100, 1000, 8000]) 
% xlabel(all_axes, 'Frequency: [Hz]'); 
% set(gcf, 'Position', [100, 100, 1200, 800]); % figure size
% % if isempty(scale) || isnan(scale)
% %     h = sgtitle(eval_params.algorithm_name);
% % else
% %     h = sgtitle( sprintf('%s (scale = %g)', eval_params.algorithm_name, scale) );
% % end
% h.FontSize  = 14;           % 字号
% h.FontName  = 'Times New Roman';  % 字体
% h.FontWeight = 'bold';      % 加粗: 'normal' / 'bold'
% h.Color     = [0 0 1];      % 颜色 (RGB)，这里是蓝色
% h.Interpreter = 'latex';  % 需要 LaTeX/TeX 语法时
% 
% 
% % BZ SPL plot
% figure;semilogx(freq_params.target_freqs,performance.BZ_SPL);title('BZ-SPL Performance');
% xlim([freq_params.target_freq_start freq_params.target_freq_end]);
% ylim([0 80]);
% 
% f_target = 1e3;
% 
% 
% % fprintf('Maximum frequency value calculated: %.2f Hz\n', f_max_val);
% 
% %% what does worst case in SFC mean ?
% 
% targe_freq = [500,1000,2000];
% for freq_num = 1:length(targe_freq)
%     figure;
%     for algorithms = 1: length(algorithms_to_test)
%         semilogx(para.temperature_vector_celsius,performance.(algorithmname).AC(:,f_idx));
% 
%     end
% 
% end