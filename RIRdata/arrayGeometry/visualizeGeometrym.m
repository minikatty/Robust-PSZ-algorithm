%% ============ 4. 可视化 ============
clc;close;
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
    view(22, 20);  % 调整观察视角
    
    % 绘制房间边界
    plotRoom(room_size);
    
    % 绘制扬声器
    hSpeaker = plot3(speakers(:,1), speakers(:,2), speakers(:,3), 'ks', ...
        'MarkerSize', 8,'MarkerEdgeColor','#000000','MarkerFaceColor','#CC4230'); 

     % 绘制虚拟源
    hVirtualSource = plot3(speakers(13,1), speakers(13,2), speakers(13,3), 's', ...
    'MarkerSize', 8, 'MarkerFaceColor', '#FFD700', 'MarkerEdgeColor', 'k');
    
    % 标注第1和第48个扬声器
    text(speakers(1,1)-0.2, speakers(1,2), speakers(1,3)+0.2, '1st', ...
        'Color', 'blue', 'FontSize', 12, 'FontName','Times New Roman', ...
        'FontWeight', 'bold');
    text(speakers(48,1)-0.3, speakers(48,2), speakers(48,3)+0.1, '48th', ...
        'Color', 'blue', 'FontSize', 12, 'FontName','Times New Roman', ...
        'FontWeight', 'bold');
    
    % 绘制暗区(DZ)
    plotCircularZone(center_DZ, r_D, [0.3 0.3 0.3]);
    hDZ_control = plot3(DZ_control(:,1), DZ_control(:,2), DZ_control(:,3), 'o', ...
        'MarkerSize', 5, 'MarkerFaceColor', [0.4 0.4 0.4], ...
        'MarkerEdgeColor', [0.4 0.4 0.4], 'LineWidth', 1);
    hDZ_monitor = plot3(DZ_monitor(:,1), DZ_monitor(:,2), DZ_monitor(:,3), 'x', ...
        'MarkerSize', 5, 'Color', 'k', 'LineWidth', 1);
    
    % 绘制亮区(BZ)
    plotCircularZone(center_BZ, r_B, [1 0.8 0]);
    hBZ_control = plot3(BZ_control(:,1), BZ_control(:,2), BZ_control(:,3), 'o', ...
        'MarkerSize', 5,'MarkerEdgeColor','#5892E8', ...
        'MarkerFaceColor','#5892E8', 'LineWidth', 1);
    hBZ_monitor = plot3(BZ_monitor(:,1), BZ_monitor(:,2), BZ_monitor(:,3), 'x', ...
        'MarkerSize', 5, 'Color', '#1ABC9C', 'LineWidth', 1);
    
    xlabel('X [m]', 'FontSize', 13,'FontName','Times New Roman','FontWeight','bold');
    ylabel('Y [m]', 'FontSize', 13,'FontName','Times New Roman','FontWeight','bold');
    zlabel('Z [m]', 'FontSize', 13,'FontName','Times New Roman','FontWeight','bold');
    % title(sprintf(['\fontname{Times New Roman}{3D} \fontname{SimSun}{房间布局' ...
    %     '&房间尺寸}: [%.1f × %.1f × %.1f] m'], room_size), ...
    %     'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('\\fontname{Times New Roman}3D \\fontname{SimSun}房间布局\n尺寸：\\fontname{Times New Roman}[%.1f×%.1f×%.1f]m', ...
    room_size), 'FontSize', 14, 'FontWeight', 'bold');
    
    xlim([0, room_size(1)]);
    ylim([0, room_size(2)]);
    zlim([0, room_size(3)]);

    % 手动指定legend（混合中英文字体）
    lgd = legend([hSpeaker, hDZ_control, hDZ_monitor, hBZ_control, hBZ_monitor, hVirtualSource], ...
        {'\fontname{Times New Roman}Speakers \fontname{SimSun}扬声器', ...
         '\fontname{Times New Roman}DZ \fontname{SimSun}控制点', ...
         '\fontname{Times New Roman}DZ \fontname{SimSun}监测点', ...
         '\fontname{Times New Roman}BZ \fontname{SimSun}控制点', ...
         '\fontname{Times New Roman}BZ \fontname{SimSun}监测点', ...
         '\fontname{Times New Roman}VS \fontname{SimSun}虚拟源'}, ...
        'Location', 'northeast', 'Interpreter', 'tex');
    lgd.NumColumns = 3;
    lgd.FontSize = 11;
    lgd.FontWeight = 'bold';
    lgd.Position = [0.1592,0.3111,0.2313,0.0806];
    grid on; axis equal;  
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
    plot(speakers(:,1), speakers(:,2), 'ks', 'MarkerSize', 10, ...
        'MarkerEdgeColor','#000000','MarkerFaceColor','#CC4230');
    text(speakers(1,1)-0.25, speakers(1,2), '1st', 'Color', 'blue', ...
        'FontSize', 12, 'FontName','Times New Roman','FontWeight', 'bold');
    text(speakers(48,1)-0.30, speakers(48,2), '48th', 'Color', 'blue', ...
        'FontSize', 12, 'FontName','Times New Roman','FontWeight', 'bold');
    text(speakers(13,1), speakers(13,2)+0.10, speakers(13,3), ...
        {'\fontname{Times New Roman}Virtual Source', '\fontname{SimSun}虚拟源'}, ...
        'Interpreter', 'tex', ...
        'Color', 'blue', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom'); 

    % 绘制虚拟源 
     plot(speakers(13,1), speakers(13,2), 's', ...
    'MarkerSize', 10, 'MarkerFaceColor', '#FFD700', 'MarkerEdgeColor', 'k');

    % 绘制暗区
    rectangle('Position', [center_DZ(1)-r_D, center_DZ(2)-r_D, 2*r_D, 2*r_D], ...
        'Curvature', [1,1], 'FaceColor', [0.3 0.3 0.3], 'FaceAlpha', 0.2, ...
        'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 2);
    plot(DZ_control(:,1), DZ_control(:,2), 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', [0.4 0.4 0.4], ...
        'MarkerEdgeColor', [0.4 0.4 0.4], ...
        'LineWidth', 0.5);
    % 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', [0.4 0.4 0.4]
    plot(DZ_monitor(:,1), DZ_monitor(:,2), 'kx', 'MarkerSize', 5, ...
         'LineWidth', 1);
    text(center_DZ(1), center_DZ(2)-r_D-0.1, 'DZ', 'FontSize', 13, ...
        'FontWeight', 'bold', 'FontName','Times New Roman', ...
        'HorizontalAlignment', 'center');
    
    % 绘制亮区
    rectangle('Position', [center_BZ(1)-r_B, center_BZ(2)-r_B, 2*r_B, 2*r_B], ...
        'Curvature', [1,1], 'FaceColor', [1 0.8 0], 'FaceAlpha', 0.2,...
        'EdgeColor', [1 0.8 0], 'LineWidth', 2);

    plot(BZ_control(:,1), BZ_control(:,2), 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', '#5892E8', 'MarkerEdgeColor', '#5892E8', 'LineWidth', 0.5);
    plot(BZ_monitor(:,1), BZ_monitor(:,2), 'x', 'MarkerSize', 5, ...
        'Color', '#1ABC9C', 'LineWidth', 1.5);
    text(center_BZ(1), center_BZ(2)-r_B-0.1, 'BZ', 'FontSize', 12, ...
        'FontWeight', 'bold', 'FontName','Times New Roman', ...
        'HorizontalAlignment', 'center');
    
    % 标注尺寸
    % 标注位置：两个圆圈下方
    y_line = center_DZ(2) - r_D - 0.18;  % 调整这个值控制高度
    plot([center_DZ(1), center_BZ(1)], [y_line, y_line], 'k-', 'LineWidth', 1);
    plot([center_DZ(1), center_DZ(1)], [y_line-0.02, y_line+0.02], 'k-', 'LineWidth', 1);
    plot([center_BZ(1), center_BZ(1)], [y_line-0.02, y_line+0.02], 'k-', 'LineWidth', 1);
    text((center_DZ(1)+center_BZ(1))/2, y_line-0.06, y_line+0.02, '1.2m', 'FontSize', 10, ...
        'HorizontalAlignment', 'center');
   
    xlabel('X [m]', 'FontSize', 13,'FontName','Times New Roman','FontWeight','bold');
    ylabel('Y [m]', 'FontSize', 13,'FontName','Times New Roman','FontWeight','bold');
    title(sprintf(['\\fontname{SimSun}俯视图 (\\fontname{Times New Roman}z = %.1f m' ...
        '\\fontname{SimSun}平面)'], center_DZ(3)), ...
        'FontSize', 14, 'FontWeight', 'bold');
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

%% ============ 绘制圆形控制区域 ===========
function plotCircularZone(center, radius, color)
    theta = linspace(0, 2*pi, 100);
    x = center(1) + radius * cos(theta);
    y = center(2) + radius * sin(theta);
    z = repmat(center(3), size(theta));
    
    fill3(x, y, z, color, 'FaceAlpha', 0.2, 'EdgeColor', color, 'LineWidth', 1);
end