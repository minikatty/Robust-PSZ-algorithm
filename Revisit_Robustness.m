clear; clc; close all;

% --- 1. path adding & data loading ---
current_mode = 'temperature';
data_dir = fullfile('data', 'SimulateRIR', current_mode);
para_file = fullfile(data_dir, 'para.mat');
load(para_file);

Temperature = para.temperature_vector_celsius; % X：51 tmperatures
freqs = para.freq_params.target_freqs;

addpath(genpath('results/temperature_AnalysisData'));
data = load('performance_matrices_20260428_0118.mat');
val_results = data.val_results;

% get uncertainty bound data
design_env_para = 22.5; % temperature
design_filename = get_data_filename(data_dir, current_mode, design_env_para);
[gamma_all, epsilon_all] = get_bound_paras_all(current_mode, design_filename);

% --- find the corresponding Index ---
target_frequency = [500, 1000, 3000];
f_indices = zeros(1, length(target_frequency));
for i = 1:length(target_frequency)
    [~, f_indices(i)] = min(abs(freqs - target_frequency(i)));
end

algorithms_to_test = {'ACC', 'PM','wcRACC', 'POTDC_RACC','RPM','RACC_PM_Subpro'};

base_colors =[
    0, 114, 189;   % ACC
    217, 83, 25;   % PM
    237, 177, 32;  % wcRACC
    126, 47, 142;  % POTDC_RACC
    119, 172, 48;  % RPM
    226, 17, 108   % RACC_PM_Subpro (Proposed)
] / 255;

line_styles = {'--', '--', '-', '-', '-', '-'}; 
markers = {'o', 's', '^', 'v', 'd', 'p'};
marker_step = 5; 
metrics_to_evaluate = {'AC'};
nominal_temp = design_env_para;

for m = 1:length(metrics_to_evaluate)
    metric_name = metrics_to_evaluate{m};       
    for f = 1:length(target_frequency)
        figure;
        hold on; grid on;        
        current_f_idx = f_indices(f);
        actual_f = freqs(current_f_idx);
        
        for alg_idx = 1:length(algorithms_to_test)
            alg_name = algorithms_to_test{alg_idx};
            y_data = val_results.(alg_name).(metric_name)(:, current_f_idx);            
            plot(Temperature, y_data, ...
                'LineStyle', line_styles{alg_idx}, ...
                'Marker', markers{alg_idx}, ...
                'Color', base_colors(alg_idx, :), ...
                'LineWidth', 1.8, ...
                'MarkerSize', 5, ...
                'MarkerIndices', 1:marker_step:length(Temperature), ... 
                'MarkerFaceColor', base_colors(alg_idx, :), ...
                'DisplayName', strrep(alg_name, '_', '-'));             
        end
        xlabel(' $||\Delta \mathbf{H}_Z||_F$ \& Temperature[$^\circ$ C]', 'FontSize', 12, 'Interpreter', 'latex');
        if f == 1
            ylabel(sprintf('%s: [dB]', metric_name), 'FontSize', 12);
        end
        if f == 2 % mark the nominal temperature
            text(nominal_temp, yl(1) + 0.95*(yl(2)-yl(1)), ' Nominal Temperature','VerticalAlignment','bottom', ...
                'FontSize', 10, 'Color', [0.3 0.3 0.3], 'FontWeight', 'bold','FontName','Times New Roman');
        end    
        set(gca, 'FontSize', 11, 'LineWidth', 1, 'Box', 'on');
        if f == 3
            legend('Location', 'best', 'FontSize', 10, 'NumColumns', length(algorithms_to_test));
        end
        ax1 = gca;
        ax1.Box = 'off';
        norm_delta_H = epsilon_all.B(current_f_idx, :); 
        xticks_bottom = xticks(ax1);        
        % 计算这些刻度位置对应的 ||Delta H||_F 数值
        % 使用 interp1 确保即使刻度没有刚好落在 51 个温度点上，也能精准插值
        top_tick_vals = interp1(Temperature, norm_delta_H, xticks_bottom, 'linear', 'extrap');
        
        % 将数值格式化为保留两位小数的字符串
        top_tick_labels = arrayfun(@(x) sprintf('%.2f', x), top_tick_vals, 'UniformOutput', false);        
        % 位置和大小与 ax1 绝对重合，但背景透明，X轴在顶部
        ax2 = axes('Position', ax1.Position, ...
                   'XAxisLocation', 'top', ...
                   'YAxisLocation', 'right', ...
                   'Color', 'none', ...      
                   'YTick',[], ...           
                   'XLim', ax1.XLim);                      
        ax2.XTick = xticks_bottom;
        ax2.XTickLabel = top_tick_labels;       
        set(ax2, 'FontSize', 11, 'LineWidth', 1);        
        linkaxes([ax1, ax2], 'x');
        xlim([min(Temperature), max(Temperature)]);
        yl = ylim;
        hold on;
        plot([nominal_temp, nominal_temp], yl, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end

