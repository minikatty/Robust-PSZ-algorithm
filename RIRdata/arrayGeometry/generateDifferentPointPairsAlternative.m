function [control_points, monitor_points] = generateDifferentPointPairsAlternative(center, radius, N)
    % 使用不同的采样策略生成两组点
    % ============ 备选方法：使用两种不同的生成模式 ============
    % 控制点：使用Vogel's螺旋（均匀分布）
    control_points_2d = generateVogelSpiral(center(1:2), radius, N);
    
    % 监测点：使用偏移的Vogel's螺旋或随机分布
    % 方法1：旋转偏移
    offset_angle = pi / 3;  % 旋转偏移角度
    monitor_points_2d = zeros(N, 2);
    for i = 1:N
        r = sqrt(sum((control_points_2d(i,:) - center(1:2)).^2));
        theta = atan2(control_points_2d(i,2) - center(2), control_points_2d(i,1) - center(1));
        theta_new = theta + offset_angle;
        
        monitor_points_2d(i, 1) = center(1) + r * cos(theta_new);
        monitor_points_2d(i, 2) = center(2) + r * sin(theta_new);
    end
    
    % 添加z坐标
    control_points = [control_points_2d, repmat(center(3), N, 1)];
    monitor_points = [monitor_points_2d, repmat(center(3), N, 1)];
end
function points_2d = generateVogelSpiral(center_2d, radius, N)
    % Vogel's 黄金角螺旋
    golden_angle = pi * (3 - sqrt(5));
    points_2d = zeros(N, 2);
    
    for i = 1:N
        theta = i * golden_angle;
        r = radius * sqrt(i / N);
        points_2d(i, 1) = center_2d(1) + r * cos(theta);
        points_2d(i, 2) = center_2d(2) + r * sin(theta);
    end
end