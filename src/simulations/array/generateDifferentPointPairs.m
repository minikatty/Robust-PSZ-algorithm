%% ============ 生成不重复的控制点和监测点对 ============
function [control_points, monitor_points] = generateDifferentPointPairs(center, radius, N)
    % 在同一圆形区域内生成两组不重复的点
    %
    % 输入:
    %   center: 圆心 [x, y, z]
    %   radius: 半径
    %   N: 每组点的数量
    %
    % 输出:
    %   control_points: 控制点 [N × 3]
    %   monitor_points: 监测点 [N × 3]
    
    % 策略：在圆内生成2N个点，然后分成两组
    % 使用改进的方法确保良好的空间分布
    
    min_distance = radius * 0.05;  % 点之间的最小距离（避免太近）
    max_attempts = 10000;
    
    all_points_2d = [];
    attempts = 0;
    
    % 生成2N个不重复的点
    while size(all_points_2d, 1) < 2*N && attempts < max_attempts
        % 使用随机采样
        r = radius * sqrt(rand());
        theta = 2 * pi * rand();
        
        x = center(1) + r * cos(theta);
        y = center(2) + r * sin(theta);
        new_point = [x, y];
        
        % 检查与已有点的距离
        if isempty(all_points_2d)
            all_points_2d = new_point;
        else
            distances = sqrt(sum((all_points_2d - new_point).^2, 2));
            if min(distances) >= min_distance
                all_points_2d = [all_points_2d; new_point];
            end
        end
        
        attempts = attempts + 1;
    end
    
    % 如果点数不够，降低距离要求
    if size(all_points_2d, 1) < 2*N
        warning('无法生成足够的分散点，降低距离要求');
        min_distance = min_distance * 0.5;
        
        while size(all_points_2d, 1) < 2*N
            r = radius * sqrt(rand());
            theta = 2 * pi * rand();
            x = center(1) + r * cos(theta);
            y = center(2) + r * sin(theta);
            new_point = [x, y];
            
            if isempty(all_points_2d)
                all_points_2d = new_point;
            else
                distances = sqrt(sum((all_points_2d - new_point).^2, 2));
                if min(distances) >= min_distance
                    all_points_2d = [all_points_2d; new_point];
                end
            end
        end
    end
    
    % 随机分成两组
    indices = randperm(size(all_points_2d, 1));
    control_indices = indices(1:N);
    monitor_indices = indices(N+1:2*N);
    
    % 添加z坐标
    control_points = [all_points_2d(control_indices, :), repmat(center(3), N, 1)];
    monitor_points = [all_points_2d(monitor_indices, :), repmat(center(3), N, 1)];
end
