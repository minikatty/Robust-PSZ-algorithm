function w = RACC_PM(H_B, DeltaH_B, H_D, DeltaH_D, p_Bd, epsilon_B, epsilon_D, para)
    %% 交替最小化迭代主循环
    max_iter = para.max_iter;   % 最大迭代次数
    tol = para.tol;       % 收敛阈值
    e_w = para.e_w;
    w_prev = para.w_prev;
    power_ratio = para.power_ratio;
    for iter = 1 : max_iter
        % ---------- Step 1: 固定误差，优化 w ---------- %
        w_curr = solve_w(H_B, DeltaH_B, H_D, DeltaH_D, p_Bd, e_w, power_ratio);
        
        % ---------- Step 2: 固定 w，优化误差 DeltaH_B、DeltaH_D ---------- %
        DeltaH_B = solve_deltaHB(H_B, w_curr, p_Bd, epsilon_B);
        DeltaH_D = solve_deltaHD(H_D, w_curr, epsilon_D);
        
        % ---------- 收敛判定 ---------- %
        if norm(w_curr - w_prev, 'fro') < tol && ...
           norm(DeltaH_B, 'fro') < tol && ...
           norm(DeltaH_D, 'fro') < tol
            fprintf('迭代 %d 次收敛！\n', iter);
            break;
        end
        w_prev = w_curr; % 更新 w
    end
    
    if iter == max_iter
        warning('达到最大迭代次数，可能未完全收敛！');
    end
    w = w_curr;
    % 3. 输出结果
    disp('最终波束赋形向量 w：');
    disp(w_curr);
    disp('最终信道误差 DeltaH_B：');
    disp(DeltaH_B);
    disp('最终信道误差 DeltaH_D：');
    disp(DeltaH_D);
end

%% 子函数 1: 固定误差时优化 w（Step 1）
function w = solve_w(H_B, DeltaH_B, H_D, DeltaH_D, p_Bd, e_w, power_ratio)
    % 构造含误差的信道矩阵（复数）
    H_tildeB = H_B + DeltaH_B;
    H_tildeD = H_D + DeltaH_D;
    
    % 拆分复数矩阵为实部和虚部（后续计算用实数运算）
    H_tildeB_real = real(H_tildeB);
    H_tildeB_imag = imag(H_tildeB);
    p_Bd_real = real(p_Bd);
    p_Bd_imag = imag(p_Bd);
    
    % 目标函数：输入实数向量 x = [w_real; w_imag]，输出实数标量
    obj_fun = @(x) norm( ...
        [p_Bd_real - (H_tildeB_real*x(1:end/2) - H_tildeB_imag*x(end/2+1:end)); ...
         p_Bd_imag - (H_tildeB_real*x(end/2+1:end) + H_tildeB_imag*x(1:end/2))] ...
    ).^2;
    
    % 约束函数：拆分 w'*R_tildeB*w 和 w'*R_tildeD*w 为实数运算
    nonlcon = @(x) power_ratio_constraint(x, H_tildeB, H_tildeD, power_ratio, e_w);
    
    % 初始猜测：实数向量 [w_real; w_imag]
    n = size(H_B, 2); % w 的维度
    w0_real = randn(n, 1); % 实部初始值（随机实数）
    w0_imag = randn(n, 1); % 虚部初始值（随机实数）
    x0 = [w0_real; w0_imag]; % 合并为实数优化变量
    
    % fmincon 优化选项
    options = optimoptions('fmincon', ...
        'Display', 'off', ...  % 关闭迭代显示
        'MaxFunctionEvaluations', 1e5, ... % 最大函数评估次数
        'MaxIterations', 1e3); % 最大迭代次数
    
    % 调用 fmincon 求解（实数优化）
    x_sol = fmincon(obj_fun, x0, [], [], [], [], [], [], nonlcon, options);
    
    % 合并实部和虚部，还原复数解
    w = x_sol(1:n) + 1i*x_sol(n+1:end);
end

% 约束函数：处理实数向量 x 的约束（拆分实虚部计算）
function [c, ceq] = power_ratio_constraint(x, H_tildeB, H_tildeD, power_ratio, e_w)
    n = length(x)/2;
    w_real = x(1:n);   % 优化变量的实部
    w_imag = x(n+1:end); % 优化变量的虚部
    
    % 约束 1: ||w||² ≤ e_w  --> w_real'*w_real + w_imag'*w_imag ≤ e_w
    c(1) = w_real'*w_real + w_imag'*w_imag - e_w;
    
    % 计算 R_tildeB = H_tildeB'*H_tildeB（复数）
    R_tildeB = H_tildeB' * H_tildeB;
    R_tildeD = H_tildeD' * H_tildeD;
    
    % 拆分 R_tildeB 和 R_tildeD 为实部和虚部
    R_tildeB_real = real(R_tildeB);
    R_tildeB_imag = imag(R_tildeB);
    R_tildeD_real = real(R_tildeD);
    R_tildeD_imag = imag(R_tildeD);
    
    % 计算 w'*R_tildeB*w（复数运算转实数计算）
    % 公式：(a+bi)'*(c+di)*(a+bi) = (a*c*a - a*d*b - b*c*b - b*d*a) + i(...)，这里取实部
    numerator_real = w_real'*R_tildeB_real*w_real - w_imag'*R_tildeB_imag*w_imag ...
                     - w_real'*R_tildeB_imag*w_imag - w_imag'*R_tildeB_real*w_real;
    % 实际更简单的方式：用复数计算再取实部
    w_complex = w_real + 1i*w_imag;
    numerator = real(w_complex' * R_tildeB * w_complex);
    denominator = real(w_complex' * R_tildeD * w_complex);
    
    % 约束 2: (w'*R_tildeB*w)/(w'*R_tildeD*w) = power_ratio
    ceq(1) = numerator - power_ratio * denominator;
    
    % 若分母可能为0，需加保护（如：denominator = max(denominator, 1e-8);）
end

%% 子函数 2: 固定 w 时优化 DeltaH_B（最坏情况误差，Step 2）
function DeltaH_B = solve_deltaHB(H_B, w, p_Bd, epsilon_B)
    e_B = p_Bd - H_B * w;  % 误差基准向量
    if norm(e_B) < 1e-9  % e_B 接近0时，误差不影响目标函数
        DeltaH_B = zeros(size(H_B), 'complex');
        return;
    end
    % 构造使误差最大的 DeltaH_B：DeltaH_B = epsilon_B * (e_B * w') / (norm(e_B)*norm(w))
    unit_eB = e_B / norm(e_B);
    unit_w = w / norm(w);
    DeltaH_B = epsilon_B * (unit_eB * unit_w');
end


%% 子函数 3: 固定 w 时优化 DeltaH_D（最坏情况误差，Step 2）
function DeltaH_D = solve_deltaHD(H_D, w, epsilon_D)
    % 构造使 ||(H_D + DeltaH_D)*w||_2 最大的 DeltaH_D（简化版，可根据需求调整）
    if norm(w) < 1e-9
        DeltaH_D = zeros(size(H_D), 'complex');
        return;
    end
    % 取随机方向（或更严谨的最优方向，这里简化为随机最大化范数）
    % 更优方式：DeltaH_D = epsilon_D * (v * w') / norm(w)，v 是 H_D*w 的方向
    v = (H_D * w) / norm(H_D * w);
    DeltaH_D = epsilon_D * (v * (w'/norm(w)));
end