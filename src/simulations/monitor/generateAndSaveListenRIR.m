function generateAndSaveListenRIR(array_file, varargin)
    % 生成RIR并保存
    %
    % 输入:
    %   array_file: 阵列信息文件（如 'array.mat'）
    %               必须包含: speakers, room_size, z_height
    %   varargin: 可选参数
    %       'grid_spacing' - 网格间距 [m]，默认 0.05
    %       'temperature' - 温度 [°C]，默认 20
    %       'beta' - 墙面反射系数，默认 0.3
    %       'fs' - 采样率 [Hz]，默认 16000
    %       'save_filename' - 保存文件名，默认 'gridRIR_data.mat'
    %       'margin' - 边界留白 [m]，默认 0.05
    
    %% 解析参数
    p = inputParser;
    addParameter(p, 'grid_spacing', 0.05, @isnumeric);
    addParameter(p, 'temperature', 20, @isnumeric);
    addParameter(p, 'beta', 0.3, @isnumeric);
    addParameter(p, 'fs', 16000, @isnumeric);
    addParameter(p, 'save_filename', 'MonitorGridRIR_data.mat', @ischar);
    addParameter(p, 'margin', 0.05, @isnumeric);
    addParameter(p,'z_height', 1.6, @isnumeric);
    parse(p, varargin{:});
    
    grid_spacing = p.Results.grid_spacing;
    temperature = p.Results.temperature;
    beta = p.Results.beta;
    fs = p.Results.fs;
    save_filename = p.Results.save_filename;
    margin = p.Results.margin;
    z_height = p.Results.z_height;
    
    fprintf('\n========== RIR 生成 ==========\n');
    
    %% 加载阵列数据
    array_data = load(array_file);
    speakers = array_data.array.s;
    room_size = array_data.array.roomSize;
    
    fprintf('阵列文件: %s\n', array_file);
    fprintf('扬声器数量: %d\n', size(speakers, 1));
    
    %% 计算声速
    c = temp2speed(temperature);
    fprintf('温度: %.1f°C, 声速: %.1fm/s\n', temperature, c);
    
    %% 生成网格
    x_vec = margin : grid_spacing : (room_size(1) - margin);
    y_vec = margin : grid_spacing : (room_size(2) - margin);
    [X, Y] = meshgrid(x_vec, y_vec);
    grid_points = [X(:), Y(:), repmat(z_height, numel(X), 1)];
    
    grid_info.x_vec = x_vec;
    grid_info.y_vec = y_vec;
    grid_info.nx = length(x_vec);
    grid_info.ny = length(y_vec);
    grid_info.X = X;
    grid_info.Y = Y;
    
    fprintf('网格: %d×%d=%d点 (间距%.1fcm)\n', ...
        length(x_vec), length(y_vec), size(grid_points, 1), grid_spacing*100);
    
    %% 计算RIR
    N_spk = size(speakers, 1);
    N_grid = size(grid_points, 1);
    total = N_spk * N_grid;
    
    fprintf('开始计算 %d 个RIR...\n', total);
    
    % 测试RIR确定长度
    truncated_time = 128; % unit: ms, include early reverberation
    N_samples = floor(truncated_time/1e3 * fs); % length of RIR
    
    % 预分配
    RIR = zeros(N_spk, N_grid, N_samples, 'single');
    
    % 批量计算
    tic;
    count = 1;
    for i = 1:N_spk
        for j = 1:N_grid
            RIR(i, j, :) = rir_generator(c, fs, grid_points(j,:), speakers(i,:), ...
                room_size, beta, N_samples);
   
            % 进度显示
            if mod(count, 3000) == 0 || count == total
                fprintf('  进度: %.1f%% (%d/%d)\n', count/total*100, count, total);
            end
            count = count + 1;
        end
    end
    
    elapsed = toc;
    fprintf('完成！耗时: %.1f分钟\n', elapsed/60);
    
    %% 保存
    save(save_filename, 'RIR', 'grid_points', 'grid_info', ...
         '-v7.3');
    
    file_info = dir(save_filename);
    fprintf('已保存: %s (%.1fMB)\n', save_filename, file_info.bytes/1024^2);
    fprintf('==============================\n\n');
end