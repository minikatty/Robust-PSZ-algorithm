function [roomArray]= generateSoundFieldGeometry()
    % 生成3D房间中的声场控制系统几何构型
    % 
    % 房间尺寸: [4, 4, 3] 米 (长×宽×高)
    % 所有扬声器和麦克风高度: z = 1.6m
    % 每个区域: 96个控制点 + 96个监测点（位置不同）
    % 
    % 输出:
    %   speakers: 扬声器坐标 [48 × 3]
    %   DZ_control: 暗区控制点 [96 × 3]
    %   BZ_control: 亮区控制点 [96 × 3]
    %   DZ_monitor: 暗区监测点 [96 × 3]
    %   BZ_monitor: 亮区监测点 [96 × 3]
    
    %% ============ 参数设置 ============
    % 房间参数
    room_size = [4, 4, 3];        % 房间尺寸 [长×宽×高] [m]
    z_height = 1.6;               % 扬声器和麦克风高度 [m]
    
    % 扬声器参数
    N_speakers = 48;              % 扬声器数量
    r_L = 1.5;                    % 扬声器阵列半径 [m]
    room_center = [room_size(1)/2, room_size(2)/2, z_height];  % 房间中心 [2, 2, 1.6]
    
    % 区域参数（在xy平面上）
    r_D = 0.2;                    % 暗区半径 [m]
    r_B = 0.2;                    % 亮区半径 [m]
    separation = 1.2;             % 两区域中心距离 [m]
    
    % DZ和BZ的中心在xy平面上，z坐标都是1.6m
    center_DZ = [room_center(1) - separation/2, room_center(2), z_height];  % [1.4, 2.0, 1.6]
    center_BZ = [room_center(1) + separation/2, room_center(2), z_height];  % [2.6, 2.0, 1.6]
    
    % 点数设置
    N_points_per_zone = 96;       % 每个区域的麦克风数(控制点和检测点一样多)
    
    %% ============ 1. 生成扬声器坐标 (3D) ============
    % 扬声器在高度z_height的平面上，围成圆形阵列
    % 第1个扬声器在最左侧（180度），clockwise排列
    theta_start = pi;
    theta_speakers = theta_start + linspace(0, -2*pi, N_speakers+1);
    theta_speakers = theta_speakers(1:end-1);
    
    speakers = zeros(N_speakers, 3);
    speakers(:, 1) = room_center(1) + r_L * cos(theta_speakers);  % x
    speakers(:, 2) = room_center(2) + r_L * sin(theta_speakers);  % y
    speakers(:, 3) = z_height;                                    % z
    
    fprintf('========== 3D房间声场几何构型 ==========\n');
    fprintf('房间尺寸: [%.1f × %.1f × %.1f] m\n', room_size);
    fprintf('工作平面高度: z = %.1f m\n\n', z_height);
    
    fprintf('扬声器:\n');
    fprintf('  数量: %d\n', N_speakers);
    fprintf('  阵列半径: %.2f m\n', r_L);
    fprintf('  第1个位置: (%.3f, %.3f, %.3f)\n', speakers(1,:));
    fprintf('  第48个位置: (%.3f, %.3f, %.3f)\n', speakers(48,:));
    
    %% ============ 2. 生成暗区(DZ)的控制点和监测点 (3D) ============
    % 控制点和监测点在xy平面上不同，但都在同一圆形区域内
    % [DZ_control, DZ_monitor] = generateDifferentPointPairs(center_DZ, r_D, N_points_per_zone);
    [DZ_control, DZ_monitor] = generateDifferentPointPairsAlternative(center_DZ, r_D, N_points_per_zone);

    fprintf('\n暗区(DZ):\n');
    fprintf('  中心位置: (%.2f, %.2f, %.2f)\n', center_DZ);
    fprintf('  半径: %.2f m\n', r_D);
    fprintf('  控制点数量: %d\n', size(DZ_control, 1));
    fprintf('  监测点数量: %d\n', size(DZ_monitor, 1));
    
    % 验证点是否重复
    min_dist_DZ = min(pdist2(DZ_control(:,1:2), DZ_monitor(:,1:2)), [], 'all');
    fprintf('  控制点与监测点最小距离: %.4f m\n', min_dist_DZ);
    
    %% ============ 3. 生成亮区(BZ)的控制点和监测点 (3D) ============
    % [BZ_control, BZ_monitor] = generateDifferentPointPairs(center_BZ, r_B, N_points_per_zone);
    [BZ_control, BZ_monitor] = generateDifferentPointPairsAlternative(center_BZ, r_B, N_points_per_zone);
    
    roomArray.s = speakers;
    roomArray.bCtrPtsPositions = BZ_control;
    roomArray.bPerPtsPositions = BZ_monitor;
    roomArray.dCtrPtsPositions = DZ_control;
    roomArray.dPerPtsPositions = DZ_monitor;
    roomArray.roomSize = room_size;
    
    fprintf('\n亮区(BZ):\n');
    fprintf('  中心位置: (%.2f, %.2f, %.2f)\n', center_BZ);
    fprintf('  半径: %.2f m\n', r_B);
    fprintf('  控制点数量: %d\n', size(BZ_control, 1));
    fprintf('  监测点数量: %d\n', size(BZ_monitor, 1));
    
    min_dist_BZ = min(pdist2(BZ_control(:,1:2), BZ_monitor(:,1:2)), [], 'all');
    fprintf('  控制点与监测点最小距离: %.4f m\n', min_dist_BZ);
    
    fprintf('\n总计:\n');
    fprintf('  扬声器: %d\n', N_speakers);
    fprintf('  DZ控制点: %d,  DZ监测点: %d\n', size(DZ_control, 1), size(DZ_monitor, 1));
    fprintf('  BZ控制点: %d,  BZ监测点: %d\n', size(BZ_control, 1), size(BZ_monitor, 1));
    fprintf('========================================\n');    
end