function [w_opt] = RPM(HB_ctrl, HD_ctrl, H_desired, para)
% optimization:
%   min  || [pB_d; 0] - [HB; sqrt(mu)*HD] * w ||_2 + epsilon_total * ||w||_2
%   s.t. ||w||_2^2 <= ew
    [~, L, num_frePoint] = size(HB_ctrl);
    scaling_factor = para.scale; % retained for backward compatibility
    has_explicit_radii = isfield(para, 'epsilon') && ...
        isfield(para.epsilon, 'B') && isfield(para.epsilon, 'D');
    mu = 1; % ACC-PM = PM: mu = 1
    pB_d = H_desired;
    wRPM = zeros(L, num_frePoint);
    % phys_cores = feature('numcores'); 
    % safe_workers = phys_cores - 2; % server: 14 local:8
    % 
    % if isempty(gcp('nocreate'))
    %     parpool('local', safe_workers); 
    % end
    for i = 1:num_frePoint % for debug
    % parfor i = 1:num_frePoint
        HB = squeeze(HB_ctrl(:, :, i));
        HD = squeeze(HD_ctrl(:, :, i));
        RB = HB'*HB;
        RD = HD'*HD;
        if has_explicit_radii
            epsilonB = local_pick_radius(para.epsilon.B, i);
            epsilonD = local_pick_radius(para.epsilon.D, i);
        else
            % Preserve the uncertainty setting used by historical scripts.
            epsilonB = 1e-2 * sqrt(norm(RB, 'fro'));
            epsilonD = 1e-2 * sqrt(norm(RD, 'fro'));
        end
        H_i = struct('B', HB_ctrl(:, :, i), 'D', HD_ctrl(:, :, i) );
        epsilon_i = struct('B', epsilonB,'D', epsilonD);%epsilonBD(i)
        pB_d_i = pB_d(:,i);
        e_w = 1e4; % squared weight-norm bound, matched to maintained RACC-PM
        wRPM(:, i) = RPM_solver(H_i, pB_d_i, epsilon_i, mu, e_w, L);
    end
    % delete(gcp('nocreate')); % end the parallel pool
    w_opt.w = wRPM;
    w_opt.scale = scaling_factor;
end

function radius = local_pick_radius(source, frequency_index)
    if isscalar(source)
        radius = source;
    elseif numel(source) >= frequency_index
        radius = source(frequency_index);
    else
        error('RPM:UncertaintyLengthMismatch', ...
            'No uncertainty radius is available for frequency index %d.', ...
            frequency_index);
    end
    validateattributes(radius, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
end

function [wRPM] = RPM_solver(H_i, pB_d, epsilon_i, mu, e_w, L)
    % HB_i (亮区传递函数), HD_i (暗区传递函数), pB_d (明区目标声场)
    % mu (权衡因子), epsilon (扰动量参数), ew (能量约束上限)
    % L (扬声器数量)
    HB_i = H_i.B;
    HD_i = H_i.D;
    epsilonD_i = epsilon_i.D;
    epsilonB_i = epsilon_i.B;

    % H_stack = [HB; sqrt(mu) * HD]:(M_B+M_D)*L
    H_stack = [HB_i; sqrt(mu) * HD_i];
    epsilon = sqrt(epsilonB_i^2 + mu * epsilonD_i^2);
    % p_stack = [pB_d; 0]:(M_B+M_D)*1
    M_D = size(HD_i, 1); 
    p_stack = [pB_d; zeros(M_D, 1)];
    
    % CVX Solver
    cvx_begin quiet
        cvx_solver sdpt3    
        variable w(L) complex   
        minimize( norm(p_stack - H_stack * w, 2) + epsilon * norm(w, 2) )
        
        subject to        
            norm(w, 2) <= sqrt(e_w);
    cvx_end

    wRPM = w;
end
