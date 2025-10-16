%% ============ 4. 可视化 ============
 load('array.mat');
 speakers = array.s;
 BZ_control = array.bCtrPtsPositions;
 BZ_monitor = array.bPerPtsPositions;
 DZ_control = array.dCtrPtsPositions;
 DZ_monitor = array.dPerPtsPositions;
 room_size = array.roomSize;
%% ============ 参数设置 ============
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

visualize3DGeometry(room_size, speakers, DZ_control, BZ_control, ...
                       DZ_monitor, BZ_monitor, center_DZ, center_BZ, r_D, r_B);

%% ============ 3D可视化函数 ============
function visualize3DGeometry(room_size, speakers, DZ_control, BZ_control, ...
                            DZ_monitor, BZ_monitor, center_DZ, center_BZ, r_D, r_B)
    
    figure('Position', [50, 50, 1400, 700]);
    
    %% 子图1：3D视图
    subplot(1, 2, 1);
    hold on; grid on; axis equal;
    view(45, 30);  % 设置视角
    
    % 绘制房间边界
    plotRoom(room_size);
    
    % 绘制扬声器
    plot3(speakers(:,1), speakers(:,2), speakers(:,3), 'ks', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'k');
    
    % 标注第1和第48个扬声器
    text(speakers(1,1)-0.3, speakers(1,2), speakers(1,3), '1st', ...
        'Color', 'blue', 'FontSize', 10, 'FontWeight', 'bold');
    text(speakers(48,1)-0.3, speakers(48,2), speakers(48,3), '48th', ...
        'Color', 'blue', 'FontSize', 10, 'FontWeight', 'bold');
    
    % 绘制暗区(DZ)
    plotCircularZone(center_DZ, r_D, [0.3 0.3 0.3]);
    plot3(DZ_control(:,1), DZ_control(:,2), DZ_control(:,3), 'o', ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.4 0.4 0.4], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    plot3(DZ_monitor(:,1), DZ_monitor(:,2), DZ_monitor(:,3), 'x', ...
        'MarkerSize', 5, 'Color', 'r', 'LineWidth', 1.5);
    
    % 绘制亮区(BZ)
    plotCircularZone(center_BZ, r_B, [1 0.8 0]);
    plot3(BZ_control(:,1), BZ_control(:,2), BZ_control(:,3), 'o', ...
        'MarkerSize', 4, 'MarkerFaceColor', [1 0.8 0], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    plot3(BZ_monitor(:,1), BZ_monitor(:,2), BZ_monitor(:,3), 'x', ...
        'MarkerSize', 5, 'Color', 'b', 'LineWidth', 1.5);
    
    xlabel('X [m]', 'FontSize', 11);
    ylabel('Y [m]', 'FontSize', 11);
    zlabel('Z [m]', 'FontSize', 11);
    title(sprintf('3D房间布局\n房间尺寸: [%.1f × %.1f × %.1f] m', room_size), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    xlim([0, room_size(1)]);
    ylim([0, room_size(2)]);
    zlim([0, room_size(3)]);
    
    legend({'', '扬声器', '', '', 'DZ控制点', 'DZ监测点', '', 'BZ控制点', 'BZ监测点'}, ...
        'Location', 'northeast', 'FontSize', 8);
    
    hold off;
    
    %% 子图2：俯视图（z = 1.6m平面）
    subplot(1, 2, 2);
    hold on; grid on; axis equal;
    
    % 绘制房间边界（俯视）
    rectangle('Position', [0, 0, room_size(1), room_size(2)], ...
        'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');
    
    % 绘制扬声器阵列圆
    theta_circle = linspace(0, 2*pi, 100);
    r_L = sqrt((speakers(1,1) - room_size(1)/2)^2 + (speakers(1,2) - room_size(2)/2)^2);
    room_center_xy = [room_size(1)/2, room_size(2)/2];
    plot(room_center_xy(1) + r_L * cos(theta_circle), ...
         room_center_xy(2) + r_L * sin(theta_circle), 'k--', 'LineWidth', 1);
    
    % 绘制扬声器
    plot(speakers(:,1), speakers(:,2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    text(speakers(1,1)-0.25, speakers(1,2), '1st', 'Color', 'blue', ...
        'FontSize', 10, 'FontWeight', 'bold');
    text(speakers(48,1)-0.25, speakers(48,2), '48th', 'Color', 'blue', ...
        'FontSize', 10, 'FontWeight', 'bold');
    
    % 绘制暗区
    rectangle('Position', [center_DZ(1)-r_D, center_DZ(2)-r_D, 2*r_D, 2*r_D], ...
        'Curvature', [1,1], 'FaceColor', [0.3 0.3 0.3 0.2], ...
        'EdgeColor', 'k', 'LineWidth', 2);
    plot(DZ_control(:,1), DZ_control(:,2), 'o', 'MarkerSize', 4, ...
        'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    plot(DZ_monitor(:,1), DZ_monitor(:,2), 'x', 'MarkerSize', 5, ...
        'Color', 'r', 'LineWidth', 1.5);
    text(center_DZ(1), center_DZ(2)-r_D-0.1, 'DZ', 'FontSize', 12, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    % 绘制亮区
    rectangle('Position', [center_BZ(1)-r_B, center_BZ(2)-r_B, 2*r_B, 2*r_B], ...
        'Curvature', [1,1], 'FaceColor', [1 1 1 0.3], ...
        'EdgeColor', 'k', 'LineWidth', 2);
    plot(BZ_control(:,1), BZ_control(:,2), 'o', 'MarkerSize', 4, ...
        'MarkerFaceColor', [1 0.8 0], 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    plot(BZ_monitor(:,1), BZ_monitor(:,2), 'x', 'MarkerSize', 5, ...
        'Color', 'b', 'LineWidth', 1.5);
    text(center_BZ(1), center_BZ(2)-r_B-0.1, 'BZ', 'FontSize', 12, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    % 标注尺寸
    plot([center_DZ(1), center_BZ(1)], [0.5, 0.5], 'k-', 'LineWidth', 1);
    plot([center_DZ(1), center_DZ(1)], [0.48, 0.52], 'k-', 'LineWidth', 1);
    plot([center_BZ(1), center_BZ(1)], [0.48, 0.52], 'k-', 'LineWidth', 1);
    text((center_DZ(1)+center_BZ(1))/2, 0.4, '1.2m', 'FontSize', 10, ...
        'HorizontalAlignment', 'center');
    
    xlabel('X [m]', 'FontSize', 11);
    ylabel('Y [m]', 'FontSize', 11);
    title(sprintf('俯视图 (z = %.1f m平面)\n○ 控制点  × 监测点', center_DZ(3)), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    xlim([0, room_size(1)]);
    ylim([0, room_size(2)]);
    
    hold off;
end

%% ============ 绘制房间边界 ============
function plotRoom(room_size)
    % 绘制3D房间的线框
    L = room_size(1);
    W = room_size(2);
    H = room_size(3);
    
    % 底面
    plot3([0 L L 0 0], [0 0 W W 0], [0 0 0 0 0], 'k-', 'LineWidth', 1.5);
    % 顶面
    plot3([0 L L 0 0], [0 0 W W 0], [H H H H H], 'k-', 'LineWidth', 1.5);
    % 垂直边
    plot3([0 0], [0 0], [0 H], 'k-', 'LineWidth', 1.5);
    plot3([L L], [0 0], [0 H], 'k-', 'LineWidth', 1.5);
    plot3([L L], [W W], [0 H], 'k-', 'LineWidth', 1.5);
    plot3([0 0], [W W], [0 H], 'k-', 'LineWidth', 1.5);
end

%% ============ 绘制圆形控制区域 ============
function plotCircularZone(center, radius, color)
    % 在3D空间中绘制圆形区域（透明圆盘）
    theta = linspace(0, 2*pi, 50);
    x_circle = center(1) + radius * cos(theta);
    y_circle = center(2) + radius * sin(theta);
    z_circle = repmat(center(3), size(theta));
    
    % 绘制圆周
    plot3(x_circle, y_circle, z_circle, 'Color', color, 'LineWidth', 2);
    
    % 绘制填充圆盘（可选）
    [X, Y] = meshgrid(linspace(center(1)-radius, center(1)+radius, 20), ...
                      linspace(center(2)-radius, center(2)+radius, 20));
    Z = ones(size(X)) * center(3);
    
    % 只保留圆内的点
    mask = (X - center(1)).^2 + (Y - center(2)).^2 <= radius^2;
    X(~mask) = NaN;
    Y(~mask) = NaN;
    Z(~mask) = NaN;
    
    surf(X, Y, Z, 'FaceColor', color, 'FaceAlpha', 0.2, ...
        'EdgeColor', 'none');
end