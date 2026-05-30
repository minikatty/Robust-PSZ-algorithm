function plot_soundfield_pressure(ATFs_monitor_data_file, eval_env_ATF_data,...
    filters_w, f_target, varargin)
% PLOT_SOUNDFIELD_PRESSURE 绘制声场声压分布的contour map
%
% 输入参数:
%   rir_data_file - 字符串, 包含RIR数据的文件路径 (e.g., 'gridRIR_data.mat')
%                   必须包含变量: RIR, grid_points, grid_info
%   filters_w          - 滤波器或权重, 支持两种格式:
%                        1. 时域滤波器: [N_speakers x filter_length] (通常为实数)
%                        2. 频域权重:   [N_speakers x 1] (通常为复数)
%   f_target      - 目标频率 (Hz)
%
%   varargin      - 可选键值对参数 ('Name', Value):
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
    addRequired(p, 'ATFs_monitor_data_file', @(x) ischar(x) || isstring(x));    
    addRequired(p, 'filters_w', @(x) isnumeric(x) || iscell(x));
    addRequired(p, 'f_target', @(x) isnumeric(x) && isscalar(x));
    
    % optional paras:
    % addParameter(p, 'fs', 16000, @isnumeric);
    % addParameter(p, 'N_fft', 2048, @isnumeric);
    % addParameter(p, 'DynamicRange', 40, @isnumeric);
    % addParameter(p, 'Colormap', 'turbo', @ischar);
    addParameter(p, 'SpeakerMasking', true, @islogical);
    addParameter(p, 'SpeakerDataFile', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'SpeakerMaskThreshold', 0.001, @isnumeric);

    addParameter(p, 'ShowZones', true, @islogical);                 % 是否显示明/暗区边界
    addParameter(p, 'BrightZoneCenter', [2.6, 2], @(x)isnumeric(x)&&numel(x)==2); % 亮区中心 [x y]
    addParameter(p, 'DarkZoneCenter',   [1.4, 2], @(x)isnumeric(x)&&numel(x)==2);% 暗区中心 [x y]
    addParameter(p, 'ZoneRadius', 0.2, @isnumeric);                 % 明/暗区半径 (m)
    addParameter(p, 'ZoneLineWidth', 1.8, @isnumeric);              % 边界线宽
    addParameter(p, 'ZoneColorBright', [0 0 0], @(x)isnumeric(x)&&numel(x)==3); % 亮区颜色
    addParameter(p, 'ZoneColorDark',   [0 0 0], @(x)isnumeric(x)&&numel(x)==3); % 暗区颜色
    addParameter(p, 'ZoneLineStyleBright', '-', @ischar);           % 亮区线型
    addParameter(p, 'ZoneLineStyleDark',   '-', @ischar);          % 暗区线型
    addParameter(p, 'ZoneLabels', false, @islogical);                % 是否在圆心处打标签
    addParameter(p, 'target_SPL', 76, @isnumeric);                % Uniform sound pressure: 76 [dB]
    
    % 执行解析
    parse(p, ATFs_monitor_data_file, filters_w, f_target, varargin{:});
    
    % 将解析结果赋给局部变量
    % fs = p.Results.fs;
    % N_fft = p.Results.N_fft;
    % dynamic_range = p.Results.DynamicRange;
    % db_colormap = p.Results.Colormap;
    SpeakerMasking = p.Results.SpeakerMasking;
    SpeakerDataFile = p.Results.SpeakerDataFile;
    ShowZones        = p.Results.ShowZones;
    BrightCenter     = p.Results.BrightZoneCenter;
    DarkCenter       = p.Results.DarkZoneCenter;
    ZoneR            = p.Results.ZoneRadius;
    ZoneLW           = p.Results.ZoneLineWidth;
    ZoneColB         = p.Results.ZoneColorBright;
    ZoneColD         = p.Results.ZoneColorDark;
    ZoneLSB          = p.Results.ZoneLineStyleBright;
    ZoneLSD          = p.Results.ZoneLineStyleDark;
    ZoneLabels       = p.Results.ZoneLabels;
    target_SPL       = p.Results.target_SPL;

    
    %% 2. 加载和准备数据
    fprintf('加载观测网格的ATFs数据: %s\n', ATFs_monitor_data_file);
    ATFs_obj = matfile(ATFs_monitor_data_file, 'Writable', false);
    monitor = ATFs_obj.monitor;
    H_all = monitor.ATF; 
    grid_points = monitor.grid_points;
    freq_params = monitor.freq_params;
    clear monitor;
    %H_all = fft(RIR, size(RIR,3), 3);  % 这里的fft参数是这样的
    % H_all_exist_in_workspace = evalin('base', "exist('H_all', 'var') ");
    % if H_all_exist_in_workspace
    %     H_all = evalin('base', 'H_all');
    % else
    %     load(grid_rir_data_file);
    % end
    % load()

    % vars_exist_in_workspace = evalin('base', "exist('RIR', 'var') " + ...
    %     "&& exist('grid_points', 'var') && exist('grid_info', 'var')");
    % 
    % if vars_exist_in_workspace
    %     % 方便画图调试
    %     fprintf('发现必需的RIR数据已在工作空间中，直接使用。\n');
    %     RIR = evalin('base', 'RIR');
    %     grid_points = evalin('base', 'grid_points');
    %     grid_info = evalin('base', 'grid_info');
    % else
    %     data = load(grid_rir_data_file);
    %     if ~isfield(data, 'RIR') || ~isfield(data, 'grid_points') || ~isfield(data, 'grid_info')
    %         error('RIR数据文件必须包含 RIR, grid_points, 和 grid_info 变量。');
    %     end
    %     RIR = data.RIR;
    %     grid_points = data.grid_points;
    %     grid_info = data.grid_info;        
    % end
    [N_positions, N_speakers, ~] = size(H_all); % 直接从文件获取 H_all 的维度
    % [N_speakers, N_positions, ~] = size(RIR);
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
    tic;
    N_fft = freq_params.nfft;
    fs = freq_params.fs;
    % --- case 1: time domain filter ---
    %这里的fft参数来自于configure_freq_parameters
    freq_axis = (0:N_fft-1) * fs / N_fft;
    [~, f_idx] = min(abs(freq_axis - f_target));
    actual_f = freq_axis(f_idx);
    freq_diff = abs(actual_f - f_target);
    % fprintf('正在计算 %.1f Hz 的声场分布...\n', actual_f);
    freq_resolution = fs / N_fft;
    fprintf('\n----------------- Monitor Point Frequency Analysis Info -----------------\n');
    fprintf('  Target Frequency:      %.2f Hz\n', f_target);
    fprintf('  Actual Frequency Used: %.2f Hz (Nearest FFT point)\n', actual_f);
    fprintf('  Difference:            %.3f Hz\n', freq_diff);
    fprintf('  FFT Resolution:        %.3f Hz (Max possible difference is ~%.3f Hz)\n', ...
        freq_resolution, freq_resolution / 2);
    fprintf('-----------------------------------------------------------\n');
    % 这里的f_idx来自于get_
    
    
    if isreal(filters_w)
        fprintf('input is frequency domain filter\n');
        W_all = fft(filters_w, N_fft, 2);
        W_f = W_all(:, f_idx);
    else
        % --- case 1: frequency domain filter ---
        fprintf('input is frequency domain filter\n');
    %     freq_params = configure_freq_parameters('fs_native',16e3, ...
    % 'target_freq_end',8e3,'rir_duration_ms',128);
        target_freqs = freq_params.target_freqs;
        [~, f_idx] = min(abs(target_freqs - f_target));
        actual_f = target_freqs(f_idx);
        freq_diff = abs(actual_f - f_target);
        % fprintf('正在计算 %.1f Hz 的声场分布...\n', actual_f);
        fprintf('\n----------------- Frequency Analysis Info -----------------\n');
        fprintf('  Target Frequency:      %.2f Hz\n', f_target);
        fprintf('  Actual Frequency Used: %.2f Hz (Nearest FFT point)\n', actual_f);
        fprintf('  Difference:            %.3f Hz\n', freq_diff);
        fprintf('-----------------------------------------------------------\n');
        W_f = filters_w(:, f_idx);
    end
    tmp_Hfs = load(eval_env_ATF_data);
    eval_ATFs_BZ = tmp_Hfs.ATF_BZ.eval;
    eval_ATF_BZ = squeeze(eval_ATFs_BZ(:, :, f_idx));
    tmp_P_complex = eval_ATF_BZ * W_f;
    p_current_RMS = sqrt(mean(abs(tmp_P_complex).^2)) + eps; % pressure density
    P_ref = 20e-6;  % ref pressure in the air
    P_target_Pa = P_ref * db2mag(target_SPL);
    alpha = P_target_Pa / p_current_RMS; % scaling factor
    W_f = W_f .* alpha;

    H_f = squeeze(H_all(:, :, f_idx));
    P_complex = H_f * W_f;    

    % P_complex = P_complex(:);
    fprintf('向量化计算完成！耗时: %.4f 秒\n', toc)
    
    %% 4. 数据后处理
    P_magnitude = abs(P_complex);
    P_phase = angle(P_complex);
    epsilon = 1e-12; 
    P_SPL = 20 * log10((P_magnitude + epsilon) / P_ref);
    
    % 重构网格
    n_rows = sqrt(length(grid_points));
    n_cols = n_rows;
    X_grid = reshape(grid_points(:, 1), n_rows, n_cols);
    X_grid = flip(X_grid);
    Y_grid = reshape(grid_points(:, 2), n_rows, n_cols);
    P_dB_grid = reshape(P_SPL, n_rows, n_cols);
    P_dB_grid(P_dB_grid > 80) = NaN; % mask the high pressure for observation
    P_phase_grid = reshape(P_phase, n_rows, n_cols);

    %% 5. 屏蔽扬声器奇点
    if SpeakerMasking
        if isempty(SpeakerDataFile)
            error("要启用 'SpeakerMasking', 必须提供 'SpeakerDataFile' 参数。");
        end
        fprintf('加载阵列布局几何数据: %s\n', SpeakerDataFile);
        array = load(SpeakerDataFile);

        fprintf('屏蔽扬声器位置奇点...\n');        
        speakers = array.roomArray.s;
        
        threshold = p.Results.SpeakerMaskThreshold;
        for i = 1:size(speakers, 1)
            distances = vecnorm(grid_points - speakers(i,:), 2, 2);
            close_points_idx = find(distances < threshold);
            if ~isempty(close_points_idx)
                % 将一维向量 P_dB 中的对应值设为 NaN
                P_SPL(close_points_idx) = NaN;
            end
        end
        % 用更新后的 P_dB 重新 reshape
        P_dB_grid = reshape(P_SPL, n_rows, n_cols);
    end
    fprintf('完成！声压级范围: %.2f ~ %.2f dB\n', ...
        min(P_SPL,[],'omitnan'), max(P_SPL,[],'omitnan'));

        %% 6. 计算区域平均声压
    if ShowZones
 
        p_ref = 2e-5; % 参考声压 (20 uPa)2e-5

        % 识别亮区内的点
        dist_from_bright = vecnorm(grid_points - [BrightCenter,1.6], 2, 2);
        bright_zone_indices = find(dist_from_bright <= ZoneR);

        % 识别暗区内的点
        dist_from_dark = vecnorm(grid_points - [DarkCenter,1.6], 2, 2);
        dark_zone_indices = find(dist_from_dark <= ZoneR);

        % 计算亮区平均平方声压并转换为dB
        P_complex_bright = P_complex(bright_zone_indices);
        mean_square_pressure_bright = sqrt(mean(abs(P_complex_bright).^2));
        avg_SPL_bright = 10 * log10(mean_square_pressure_bright / (p_ref^2));

        % 计算暗区平均平方声压并转换为dB
        P_complex_dark = P_complex(dark_zone_indices);
        mean_square_pressure_dark = mean(abs(P_complex_dark).^2);
        avg_SPL_dark = 10 * log10(mean_square_pressure_dark / (p_ref^2));

        % 计算声学对比度
        acoustic_contrast = avg_SPL_bright - avg_SPL_dark;

        % 在命令窗口显示结果
        fprintf('\n----------------- 区域声压分析结果 (@ %.1f Hz) -----------------\n', actual_f);
        fprintf('亮区 (Bright Zone) @ [%.2f, %.2f], r=%.2f:\n', BrightCenter(1), BrightCenter(2), ZoneR);
        fprintf('  - 观测点数量: %d\n', length(bright_zone_indices));
        fprintf('  - 平均声压级: %.2f dB\n', avg_SPL_bright);
        fprintf('\n暗区 (Dark Zone) @ [%.2f, %.2f], r=%.2f:\n', DarkCenter(1), DarkCenter(2), ZoneR);
        fprintf('  - 观测点数量: %d\n', length(dark_zone_indices));
        fprintf('  - 平均声压级: %.2f dB\n', avg_SPL_dark);
        fprintf('\n总体性能:\n');
        fprintf('  - 声学对比度: %.2f dB\n', acoustic_contrast);
        fprintf('--------------------------------------------------------------------\n');
    end

    %% 6. 绘图
    % 声压级(dB)分布
    figure;
    contourf(X_grid, Y_grid, P_dB_grid, 20, 'LineStyle', 'none');
    axis equal tight; grid on;
    xlabel('X (m)', 'FontSize', 12); ylabel('Y (m)', 'FontSize', 12);
    xticks(0:1:4);yticks(0:1:4);
    % title(sprintf('声压级(dB)分布 @ %.1f Hz', actual_f), 'FontSize', 14);
    % colormap(turbo);
   
    colormap(brewermap([], '-RdBu')); %inferno
    c = colorbar; c.Label.String = 'SPL (dB)'; c.Label.FontSize = 11;
    % max_dB = max(P_dB_grid(:),[],'omitnan');
    % min_dB = min(P_dB_grid(:),[],'omitnan');
    max_dB = 80;
    min_dB = 50;
    clim([min_dB,max_dB]);
    if ShowZones
        hold on;
        draw_circle(gca, BrightCenter, ZoneR, 'Color', ZoneColB, 'LineWidth', ZoneLW, 'LineStyle', ZoneLSB);
        draw_circle(gca, DarkCenter,   ZoneR, 'Color', ZoneColD, 'LineWidth', ZoneLW, 'LineStyle', ZoneLSD);
        if ZoneLabels
            text(BrightCenter(1), BrightCenter(2), 'Bright', 'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontWeight','bold','Color',ZoneColB);
            text(DarkCenter(1),   DarkCenter(2),   'Dark',   'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontWeight','bold','Color',ZoneColD);
        end
        hold off;
    end

    grid off;

    % 相位分布
    figure;
    contourf(X_grid, Y_grid, P_phase_grid, 20, 'LineStyle', 'none');
    if ShowZones
        hold on;
        draw_circle(gca, BrightCenter, ZoneR, 'Color', ZoneColB, 'LineWidth', ZoneLW, 'LineStyle', ZoneLSB);
        draw_circle(gca, DarkCenter,   ZoneR, 'Color', ZoneColD, 'LineWidth', ZoneLW, 'LineStyle', ZoneLSD);
        if ZoneLabels
            text(BrightCenter(1), BrightCenter(2), 'BZ', 'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontWeight','bold','Color',ZoneColB);
            text(DarkCenter(1),   DarkCenter(2),   'DZ',   'HorizontalAlignment','middle','FontWeight','bold','Color',ZoneColD);
        end
        hold off;
    end

    axis equal tight; grid on;
    xlim([1, 3]);
    ylim([1.5, 2.5]);
    xlabel('X (m)', 'FontSize', 12); ylabel('Y (m)', 'FontSize', 12);
    xticks(0:1:4);yticks(0:1:4);
    % title(sprintf('声压相位分布 @ %.1f Hz (rad)', actual_f), 'FontSize', 14);
    colormap(brewermap([], '-RdBu')); 
    c = colorbar; c.Label.String = '相位 (rad)'; c.Label.FontSize = 11;
    grid off;
    
end

function draw_circle(ax, center, radius, varargin)
% 以参数方程绘制圆，不依赖 toolboxes
    theta = linspace(0, 2*pi, 361);
    xx = center(1) + radius * cos(theta);
    yy = center(2) + radius * sin(theta);
    plot(ax, xx, yy, varargin{:});
end
