clear;close all;clc;
% 生成3D房间中的声场控制系统几何构型
% 
% 房间尺寸: [4, 4, 3] 米 (长×宽×高)
% 所有扬声器和麦克风高度: z = 1.6m
% 每个区域: 96个控制点 + 96个监测点（位置不同）
% Generates a dense grid of RIRs for sound field visualization.
% This is a self-contained function that intelligently manages its own parallel pool.
% It generates a large .mat file containing RIRs from all speakers to all 
% points on a 3D grid, suitable for creating SPL contour maps.
%
% 'array_file'    - string: Path to the array geometry file.
% 'output_file'   - string: Path to save the output .mat file.
% 'grid_spacing'  - scalar: Spacing for the monitor grid (m).
% 'temperature'   - scalar: Room temperature in Celsius.
% 'beta'          - scalar: Wall reflection coefficient.
% 
% 输出:
%   speakers: 扬声器坐标 [48 × 3]
%   DZ_control: 暗区控制点 [96 × 3]
%   BZ_control: 亮区控制点 [96 × 3]
%   DZ_monitor: 暗区监测点 [96 × 3]
%   BZ_monitor: 亮区监测点 [96 × 3]
%   SPL_monitor: 画声压图得监测点

addpath(genpath('src/.'));
addpath(genpath('data/.'));
addpath(genpath('lib/.'));
%% ============ 1. 基本参数设置 ============
% 房间参数
room_size = [4, 4, 3];        % 房间尺寸 [长×宽×高] [m]
z_height = 1.6;               % 扬声器和麦克风高度 [m]

% 标准环境参数
temperature = 20;             % 摄氏度
beta = 0.3;                   % 混响参数非墙体反射系数
% no noise，no position mismatch, no temperature fluctuation
    
% 扬声器参数
N_speakers = 48;              % 扬声器数量
r_L = 1.5;                    % 扬声器阵列半径 [m]
room_center = [room_size(1)/2, room_size(2)/2, z_height];  % 房间中心 [2, 2, 1.6]

% 区域参数（在xy平面上）
r_D = 0.2;                    % 暗区半径 [m]
r_B = 0.2;                    % 明区半径 [m]
separation = 1.2;             % 两区域中心距离 [m]

% DZ和BZ的中心在xy平面上，z坐标都是1.6m
center_DZ = [room_center(1) - separation/2, room_center(2), z_height];  % [1.4, 2.0, 1.6]
center_BZ = [room_center(1) + separation/2, room_center(2), z_height];  % [2.6, 2.0, 1.6]
    
% 点数设置
N_points_per_zone = 96;       % 每个区域的麦克风数(控制点和检测点一样多)
    
%% ============ 2. 生成扬声器坐标 (3D) ============
% 扬声器在高度z_height的平面上，围成圆形阵列
% 第1个扬声器在最左侧（180度），clockwise排列    
theta_start = pi;
theta_speakers = theta_start + linspace(0, -2*pi, N_speakers+1);
theta_speakers = theta_speakers(1:end-1);
    
speakers = zeros(N_speakers, 3);
speakers(:, 1) = room_center(1) + r_L * cos(theta_speakers);  % x
speakers(:, 2) = room_center(2) + r_L * sin(theta_speakers);  % y
speakers(:, 3) = z_height;                                    % z
    
% fprintf('========== 3D房间声场几何构型 ==========\n');
% fprintf('房间尺寸: [%.1f × %.1f × %.1f] m\n', room_size);
% fprintf('工作平面高度: z = %.1f m\n\n', z_height);
% fprintf('扬声器:\n');
% fprintf('  数量: %d\n', N_speakers);
% fprintf('  阵列半径: %.2f m\n', r_L);
% fprintf('  第1个位置: (%.3f, %.3f, %.3f)\n', speakers(1,:));
% fprintf('  第48个位置: (%.3f, %.3f, %.3f)\n', speakers(48,:));

    
%% ============ 3. 生成暗区(DZ)的控制点和监测点 (3D) ============
% 验证点是否重复或太近
[DZ_control, DZ_monitor] = generateDifferentPointPairs(center_DZ, r_D, N_points_per_zone);
min_dist_DZ = min(pdist2(DZ_control(:,1:2), DZ_monitor(:,1:2)), [], 'all');
% fprintf('  暗区控制点与监测点最小距离: %.4f m\n', min_dist_DZ);
[BZ_control, BZ_monitor] = generateDifferentPointPairs(center_BZ, r_B, N_points_per_zone);
min_dist_BZ = min(pdist2(BZ_control(:,1:2), BZ_monitor(:,1:2)), [], 'all');
% fprintf('  明区控制点与监测点最小距离: %.4f m\n', min_dist_BZ);


%% ============ 4.生成整个房间xy平面的监测点 (3D) ============
grid_spacing = 0.02; % 
x_vec = 0 : grid_spacing : room_size(1);
y_vec = 0 : grid_spacing : room_size(2);
[X, Y] = meshgrid(x_vec, y_vec);
grid_points = [X(:), Y(:), repmat(z_height, numel(X), 1)];

% % lsk sources and ctrl&eval&monitor mics display
% figure;scatter(grid_points(:,1), grid_points(:,2), 10); hold on;
% scatter(BZ_control(:,1),BZ_control(:,2),'rx'); hold on; % red
% scatter(BZ_monitor(:,1),BZ_monitor(:,2),'bs'); hold on; % blue
% scatter(DZ_control(:,1),DZ_control(:,2),'gx'); hold on; % green
% scatter(DZ_monitor(:,1),DZ_monitor(:,2),'ys'); hold on; % yellow
% scatter(speakers(:,1),speakers(:,2),'kd'); hold on;
% axis equal;
% xlim([0 room_size(1)]);ylim([0 room_size(2)]);

%% ============ 5.保存BZ&DZ的控制点和监测点 ============
roomArray.s = speakers;
roomArray.BZ_ctrl = BZ_control;
roomArray.BZ_eval = BZ_monitor;
roomArray.DZ_ctrl = DZ_control;
roomArray.DZ_eval = DZ_monitor;
roomArray.roomSize = room_size;
roomArray.SPL_Monitor = grid_points; % 注意画声压时的位置坐标对应

% % basic params display
% fprintf('\n明区(BZ):\n');
% fprintf('  中心位置: (%.2f, %.2f, %.2f)\n', center_BZ);
% fprintf('  半径: %.2f m\n', r_B);
% fprintf('  控制点数量: %d\n', size(BZ_control, 1));
% fprintf('  监测点数量: %d\n', size(BZ_monitor, 1));
% 
% fprintf('\n总计:\n');
% fprintf('  扬声器: %d\n', N_speakers);
% fprintf('  DZ控制点: %d,  DZ监测点: %d\n', size(DZ_control, 1), size(DZ_monitor, 1));
% fprintf('  BZ控制点: %d,  BZ监测点: %d\n', size(BZ_control, 1), size(BZ_monitor, 1));
% fprintf('========================================\n');    
save("data/arrayGeometry/array_layout.mat","roomArray");

%%  6.生成画声压图需要的数据

% 用温度计算声速
c = temp2speed(temperature);
% fprintf('温度: %.1f°C, 声速: %.1fm/s\n', temperature, c);
fprintf('\n========== Monitor_RIR 生成 ==========\n'); 
N_grid = size(grid_points, 1);

% RIR and ATFs paras settings
freq_params = configure_freq_parameters();
N_samples = freq_params.rir_len;
fs = freq_params.fs;
target_freqs = freq_params.target_freqs;
RIR = zeros(N_grid, N_speakers, N_samples, 'single');
N_freq = length(target_freqs);
ATF = zeros(N_grid, N_speakers, N_freq, 'single') + 1i;

% 开启并行池前的相关参数检查
phys_cores = feature('numcores'); 
safe_workers = max(1, phys_cores - 2); % server: 14 local:8

if isempty(gcp('nocreate'))
    parpool('local', safe_workers); 
end

% 消除阶截断带来的Gibbs效应，最后尾部加一个fade_out窗口
fade_time = 8; % 最后 8ms 加窗
fade_samples = round(fade_time / 1000 * fs); % 窗长点数
w_tail = hann(2 * fade_samples);
w_tail = w_tail(fade_samples+1 : end);

% % for debug 
% grid_points = grid_points(1:5,:);
% N_grid = size(grid_points, 1);
% RIR = zeros(N_grid, N_speakers, N_samples, 'single');
% ATF = zeros(N_grid, N_speakers, N_freq, 'single') + 1i;

% 生成monitor点位的RIR并计算控制频点的ATFs
parfor i = 1:N_speakers % parfor
    RIR_tmp = rir_generator(c, fs, grid_points, speakers(i,:), ...
        room_size, beta, N_samples);
    RIR_tmp(:,end-fade_samples+1:end) = RIR_tmp(:, end-fade_samples+1:end) .* w_tail';
    RIR_tmp = single(RIR_tmp);
    RIR_tmp = reshape(RIR_tmp, N_grid, 1, N_samples);
    RIR(:, i, :) = RIR_tmp;
    ATF_tmp = compute_atf(RIR_tmp, freq_params);
    ATF_tmp = reshape(ATF_tmp, 1, N_grid, N_freq);
    ATF(:, i, :) = ATF_tmp;
end

% 保存
monitor.grid_points = grid_points;   % 声压与该坐标对是一一对应的
monitor.ATF = ATF;
monitor.freq_params = freq_params;
save_path = "data/MonitorGrid/Monitor_IR.mat";
save(save_path, 'RIR', '-v7.3');
save_path = "data/MonitorGrid/Monitor_ATF.mat";
save(save_path, 'monitor', '-v7.3');
clear RIR monitor; %释放内存

%%  7.生成画稳健算法需要的所有数据
% --- Modes to Generate ---
modes_to_generate = {'temperature'};
% 可调节case：'snr', 'temperature', 'position'
% 'temperature'
rng(2025); % For global reproducibility

params = freq_params;
params.beta = beta;
params.w_tail = w_tail;
params.fade_samples = fade_samples;

% --- Loop through each generation mode ---
for i = 1:length(modes_to_generate)
    current_mode = modes_to_generate{i};
    generate_rir_ATF_database(current_mode, params);
end

delete(gcp('nocreate'));


