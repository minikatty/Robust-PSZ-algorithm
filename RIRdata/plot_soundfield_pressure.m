function plot_soundfield_pressure(grid_rir_data_file, ...
    filters_w, f_target, varargin)
% PLOT_SOUNDFIELD_PRESSURE 绘制声场声压分布的contour map
%
% 输入参数:
%   rir_data_file - 字符串, 包含RIR数据的文件路径 (e.g., 'gridRIR_data.mat')
%                   必须包含变量: RIR, grid_points, grid_info
%   filters_w     - [N_speakers x filter_length] 或 cell{N_speakers x 1}
%   f_target      - 目标频率 (Hz)
%
% 可选键值对参数 ('Name', Value):
%   'fs'                 - 采样率 (Hz). 默认: 16000
%   'N_fft'              - FFT点数. 默认: 2048
%   'DynamicRange'       - dB图的动态范围 (dB). 默认: 40
%   'Colormap'           - dB图的色图. 默认: 'parula'
%   'SpeakerMasking'     - 是否屏蔽扬声器位置的奇点. 默认: true
%   'SpeakerDataFile'    - 扬声器坐标文件路径 (e.g., 'array.mat'),
%                          当 'SpeakerMasking' 为 true 时必需.
%   'SpeakerMaskThreshold' - 屏蔽扬声器位置的半径 (m). 默认: 0.05
%
% example:
%   % 基础调用
%   plot_soundfield_pressure('gridRIR_data.mat', filters_w, 1000);
%
%   % 自定义参数调用
%   plot_soundfield_pressure('gridRIR_data.mat', filters_w, 1000, ...
%       'DynamicRange', 30, 'Colormap', 'viridis', 'SpeakerMasking', true, ...
%       'SpeakerDataFile', 'array.mat');

    %% 1. 解析输入参数
    p = inputParser;
    
    % 必需参数
    addRequired(p, 'grid_rir_data_file', @(x) ischar(x) || isstring(x));    
    addRequired(p, 'filters_w', @(x) isnumeric(x) || iscell(x));
    addRequired(p, 'f_target', @(x) isnumeric(x) && isscalar(x));
    
    % 可选参数
    addParameter(p, 'fs', 16000, @isnumeric);
    addParameter(p, 'N_fft', 2048, @isnumeric);
    % addParameter(p, 'DynamicRange', 40, @isnumeric);
    % addParameter(p, 'Colormap', 'turbo', @ischar);
    addParameter(p, 'SpeakerMasking', true, @islogical);
    addParameter(p, 'SpeakerDataFile', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'SpeakerMaskThreshold', 0.001, @isnumeric);
    
    % 执行解析
    parse(p, grid_rir_data_file, filters_w, f_target, varargin{:});
    
    % 将解析结果赋给局部变量
    fs = p.Results.fs;
    N_fft = p.Results.N_fft;
    % dynamic_range = p.Results.DynamicRange;
    % db_colormap = p.Results.Colormap;
    SpeakerMasking = p.Results.SpeakerMasking;
    SpeakerDataFile = p.Results.SpeakerDataFile;
    
    %% 2. 加载和准备数据
    fprintf('加载观测网格的RIR数据: %s\n', grid_rir_data_file);
    data = load(grid_rir_data_file);
    if ~isfield(data, 'RIR') || ~isfield(data, 'grid_points') || ~isfield(data, 'grid_info')
        error('RIR数据文件必须包含 RIR, grid_points, 和 grid_info 变量。');
    end
    RIR = data.RIR;
    grid_points = data.grid_points;
    grid_info = data.grid_info;
    
    [N_speakers, N_positions, ~] = size(RIR);
    fprintf('数据信息: %d 个扬声器, %d 个观测点\n', N_speakers, N_positions);
    
    % 处理滤波器格式
    if iscell(filters_w)
        filter_length = length(filters_w{1});
        filters_w_mat = zeros(N_speakers, filter_length);
        for i = 1:N_speakers
            filters_w_mat(i, :) = filters_w{i};
        end
        filters_w = filters_w_mat;
    end
    
    %% 3. 核心计算
    freq_axis = (0:N_fft-1) * fs / N_fft;
    [~, f_idx] = min(abs(freq_axis - f_target));
    actual_f = freq_axis(f_idx);
    freq_diff = abs(actual_f - f_target);
    % fprintf('正在计算 %.1f Hz 的声场分布...\n', actual_f);
    freq_resolution = fs / N_fft;
    fprintf('\n----------------- Frequency Analysis Info -----------------\n');
    fprintf('  Target Frequency:      %.2f Hz\n', f_target);
    fprintf('  Actual Frequency Used: %.2f Hz (Nearest FFT point)\n', actual_f);
    fprintf('  Difference:            %.3f Hz\n', freq_diff);
    fprintf('  FFT Resolution:        %.3f Hz (Max possible difference is ~%.3f Hz)\n', freq_resolution, freq_resolution / 2);
    fprintf('-----------------------------------------------------------\n');
    tic;
    H_all = fft(RIR, N_fft, 3);
    W_all = fft(filters_w, N_fft, 2);
    H_f = squeeze(H_all(:, :, f_idx));
    W_f = W_all(:, f_idx);
    P_complex = sum(W_f .* H_f, 1);
    P_complex = P_complex(:);
    fprintf('向量化计算完成！耗时: %.4f 秒\n', toc)
    
    %% 4. 数据后处理
    P_magnitude = abs(P_complex);
    P_phase = angle(P_complex);
    epsilon = 1e-12; 
    P_dB = 20 * log10(P_magnitude + epsilon);
    
    % 重构网格
    n_rows = grid_info.ny;
    n_cols = grid_info.nx;
    X_grid = reshape(grid_points(:, 1), n_rows, n_cols);
    Y_grid = reshape(grid_points(:, 2), n_rows, n_cols);
    P_dB_grid = reshape(P_dB, n_rows, n_cols);
    P_phase_grid = reshape(P_phase, n_rows, n_cols);
    
    %% 5. 屏蔽扬声器奇点
    if SpeakerMasking
        if isempty(SpeakerDataFile)
            error("要启用 'SpeakerMasking', 必须提供 'SpeakerDataFile' 参数。");
        end
        fprintf('加载阵列布局几何数据: %s\n', SpeakerDataFile);
        array = load(SpeakerDataFile);

        fprintf('屏蔽扬声器位置奇点...\n');        
        speakers = array.array.s;
        
        threshold = p.Results.SpeakerMaskThreshold;
        for i = 1:size(speakers, 1)
            distances = vecnorm(grid_points - speakers(i,:), 2, 2);
            close_points_idx = find(distances < threshold);
            if ~isempty(close_points_idx)
                % 将一维向量 P_dB 中的对应值设为 NaN
                P_dB(close_points_idx) = NaN;
            end
        end
        % 用更新后的 P_dB 重新 reshape
        P_dB_grid = reshape(P_dB, n_rows, n_cols);
    end
    fprintf('完成！声压级范围: %.2f ~ %.2f dB\n', ...
        min(P_dB,[],'omitnan'), max(P_dB,[],'omitnan'));

    %% 6. 绘图
    figure('Position', [100 100 1200 500]);
    % 声压级(dB)分布
    subplot(1,2,1);
    contourf(X_grid, Y_grid, P_dB_grid, 20, 'LineStyle', 'none');
    axis equal tight; grid on;
    xlabel('X (m)', 'FontSize', 12); ylabel('Y (m)', 'FontSize', 12);
    xticks(0:1:4);yticks(0:1:4);
    title(sprintf('声压级(dB)分布 @ %.1f Hz', actual_f), 'FontSize', 14);
    colormap(subplot(1,2,1), brewermap([], 'RdBu')); 
    c = colorbar; c.Label.String = '声压级 (dB)'; c.Label.FontSize = 11;
    max_dB = max(P_dB_grid(:),[],'omitnan');
    min_dB = min(P_dB_grid(:),[],'omitnan');
    clim([min_dB,max_dB]);
    grid off;

    % 相位分布
    subplot(1,2,2);
    contourf(X_grid, Y_grid, P_phase_grid, 20, 'LineStyle', 'none');
    axis equal tight; grid on;
    xlabel('X (m)', 'FontSize', 12); ylabel('Y (m)', 'FontSize', 12);
    xticks(0:1:4);yticks(0:1:4);
    title(sprintf('声压相位分布 @ %.1f Hz (rad)', actual_f), 'FontSize', 14);
    colormap(subplot(1,2,2), brewermap([], 'RdBu')); 
    c = colorbar; c.Label.String = '相位 (rad)'; c.Label.FontSize = 11;
    grid off;
    
end