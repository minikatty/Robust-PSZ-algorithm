function plot_soundfield_pressure_V2(grid_rir_data_file, ...
    filters_w, f_target, varargin)
% PLOT_SOUNDFIELD_PRESSURE 绘制声场声压分布的contour map
% ... (函数帮助文档保持不变) ...

    %% 1. 解析输入参数
    p = inputParser;
    
    % 必需参数
    addRequired(p, 'grid_rir_data_file', @(x) ischar(x) || isstring(x));    
    addRequired(p, 'filters_w', @(x) isnumeric(x) || iscell(x));
    addRequired(p, 'f_target', @(x) isnumeric(x) && isscalar(x));
    
    % 可选参数
    addParameter(p, 'fs', 16000, @isnumeric);
    addParameter(p, 'N_fft', 2048, @isnumeric);
    addParameter(p, 'SpeakerMasking', true, @islogical);
    addParameter(p, 'SpeakerDataFile', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'SpeakerMaskThreshold', 0.001, @isnumeric);
    addParameter(p, 'ShowZones', true, @islogical);
    addParameter(p, 'BrightZoneCenter', [2.6, 2], @(x)isnumeric(x)&&numel(x)==2);
    addParameter(p, 'DarkZoneCenter',   [1.4, 2], @(x)isnumeric(x)&&numel(x)==2);
    addParameter(p, 'ZoneRadius', 0.2, @isnumeric);
    addParameter(p, 'ZoneLineWidth', 1.8, @isnumeric);
    addParameter(p, 'ZoneColorBright', [0 0 0], @(x)isnumeric(x)&&numel(x)==3);
    addParameter(p, 'ZoneColorDark',   [0 0 0], @(x)isnumeric(x)&&numel(x)==3);
    addParameter(p, 'ZoneLineStyleBright', '-', @ischar);
    addParameter(p, 'ZoneLineStyleDark',   '-', @ischar);
    addParameter(p, 'ZoneLabels', true, @islogical);
    
    % 执行解析
    parse(p, grid_rir_data_file, filters_w, f_target, varargin{:});
    
    % 将解析结果赋给局部变量
    fs = p.Results.fs;
    N_fft = p.Results.N_fft;
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
    tic;
    H_all = fft(RIR, N_fft, 3);
    
    freq_axis = (0:N_fft-1) * fs / N_fft;
    [~, f_idx] = min(abs(freq_axis - f_target));
    actual_f = freq_axis(f_idx);
    
    if isreal(filters_w)
        fprintf('输入为时域滤波器...\n');
        W_all = fft(filters_w, N_fft, 2);
        W_f = W_all(:, f_idx);
    else
        fprintf('输入为频域权重...\n');
        % 假设频域权重 filters_w 是一个 [N_speakers x N_freqs] 的矩阵
        % 我们需要从中选择对应目标频率的列
        W_f = filters_w(:, f_idx); 
    end

    H_f = squeeze(H_all(:, :, f_idx));
    
    % H_f 维度: [N_speakers x N_positions]
    % W_f 维度: [N_speakers x 1]
    % P_complex 维度: [N_positions x 1]
    P_complex = H_f' * W_f;
    
    fprintf('向量化计算完成！耗时: %.4f 秒\n', toc)
    
    %% 4. 数据后处理
    P_magnitude = abs(P_complex);
    P_phase = angle(P_complex);
    epsilon = 1e-12; 
    P_dB = 20 * log10(P_magnitude + epsilon);
    
    % 重构网格
    n_rows = grid_info.ny; % 行数对应 Y
    n_cols = grid_info.nx; % 列数对应 X
    
    % X_grid 和 Y_grid 会恢复出 meshgrid 的原始结构，维度为 [ny x nx]
    X_grid = reshape(grid_points(:, 1), n_rows, n_cols);
    Y_grid = reshape(grid_points(:, 2), n_rows, n_cols);
    
    % P_dB_grid 的维度为 [ny x nx]，其行对应Y轴，列对应X轴
    P_dB_grid = reshape(P_dB, n_rows, n_cols);
    P_phase_grid = reshape(P_phase, n_rows, n_cols);
    
    %% 5. 屏蔽扬声器奇点
    if SpeakerMasking
        if isempty(SpeakerDataFile)
            error("要启用 'SpeakerMasking', 必须提供 'SpeakerDataFile' 参数。");
        end
        array = load(SpeakerDataFile);
        speakers = array.array.s;
        
        threshold = p.Results.SpeakerMaskThreshold;
        for i = 1:size(speakers, 1)
            distances = vecnorm(grid_points - speakers(i,:), 2, 2);
            close_points_idx = find(distances < threshold);
            if ~isempty(close_points_idx)
                P_dB(close_points_idx) = NaN;
            end
        end
        P_dB_grid = reshape(P_dB, n_rows, n_cols);
    end
    fprintf('完成！声压级范围: %.2f ~ %.2f dB\n', min(P_dB,[],'omitnan'), max(P_dB,[],'omitnan'));

    %% 6. 计算区域平均声压
    if ShowZones
        p_ref = 2e-5; % 标准参考声压 (20 uPa)

        % 识别亮区内的点 (在3D空间中计算距离)
        dist_from_bright = vecnorm(grid_points - [BrightCenter, grid_points(1,3)], 2, 2);
        bright_zone_indices = find(dist_from_bright <= ZoneR);

        % 识别暗区内的点 (在3D空间中计算距离)
        dist_from_dark = vecnorm(grid_points - [DarkCenter, grid_points(1,3)], 2, 2);
        dark_zone_indices = find(dist_from_dark <= ZoneR);

        % 计算亮区平均平方声压并转换为dB
        P_complex_bright = P_complex(bright_zone_indices);
        mean_square_pressure_bright = mean(abs(P_complex_bright).^2);
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

    %% 7. 绘图
    figure('Position', [100 100 1200 500]);
    % 声压级(dB)分布
    subplot(1,2,1);
    
    % --- 关键修正：对数据矩阵 P_dB_grid 进行转置(')以匹配X_grid和Y_grid的坐标 ---
    contourf(X_grid, Y_grid, P_dB_grid', 20, 'LineStyle', 'none');
    
    axis equal tight; grid on;
    xlabel('X (m)', 'FontSize', 12); ylabel('Y (m)', 'FontSize', 12);
    xticks(0:1:4);yticks(0:1:4);
    title(sprintf('声压级(dB)分布 @ %.1f Hz', actual_f), 'FontSize', 14);
    colormap(subplot(1,2,1), brewermap([], '-RdBu')); 
    c = colorbar; c.Label.String = '声压级 (dB)'; c.Label.FontSize = 11;
    max_dB = max(P_dB_grid(:),[],'omitnan');
    min_dB = min(P_dB_grid(:),[],'omitnan');
    clim([min_dB,max_dB]);
    if ShowZones
        hold on;
        draw_circle(gca, BrightCenter, ZoneR, 'Color', ZoneColB, 'LineWidth', ZoneLW, 'LineStyle', ZoneLSB);
        draw_circle(gca, DarkCenter,   ZoneR, 'Color', ZoneColD, 'LineWidth', ZoneLW, 'LineStyle', ZoneLSD);
        if ZoneLabels
            text(BrightCenter(1), BrightCenter(2), 'Bright', 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontWeight','bold','Color',ZoneColB);
            text(DarkCenter(1),   DarkCenter(2),   'Dark',   'HorizontalAlignment','center', 'VerticalAlignment','middle','FontWeight','bold','Color',ZoneColD);
        end
        hold off;
    end
    grid off;

    % 相位分布
    subplot(1,2,2);
    
    % --- 关键修正：对数据矩阵 P_phase_grid 进行转置(')以匹配X_grid和Y_grid的坐标 ---
    contourf(X_grid, Y_grid, P_phase_grid', 20, 'LineStyle', 'none');
    
    if ShowZones
        hold on;
        draw_circle(gca, BrightCenter, ZoneR, 'Color', ZoneColB, 'LineWidth', ZoneLW, 'LineStyle', ZoneLSB);
        draw_circle(gca, DarkCenter,   ZoneR, 'Color', ZoneColD, 'LineWidth', ZoneLW, 'LineStyle', ZoneLSD);
        if ZoneLabels
            text(BrightCenter(1), BrightCenter(2), 'BZ', 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontWeight','bold','Color',ZoneColB);
            text(DarkCenter(1),   DarkCenter(2),   'DZ',   'HorizontalAlignment','center', 'VerticalAlignment','middle','FontWeight','bold','Color',ZoneColD);
        end
        hold off;
    end
    axis equal tight; grid on;
    xlabel('X (m)', 'FontSize', 12); ylabel('Y (m)', 'FontSize', 12);
    xticks(0:1:4);yticks(0:1:4);
    title(sprintf('声压相位分布 @ %.1f Hz (rad)', actual_f), 'FontSize', 14);
    colormap(subplot(1,2,2), brewermap([], '-RdBu')); 
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