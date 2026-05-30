function [w_opt] = POTDC_RACC(HB_ctrl, HD_ctrl, para)

    [~, L, num_frePoint] = size(HB_ctrl);
    scaling_factor = sqrt(para.scale); % control the bound of uncertainty
    % eta = scaling_factor.* para.epsilon.B;
    % gammaD = scaling_factor.* para.gamma.D;
    wPRACC = zeros(L, num_frePoint);
    iterations = 20;
    % phys_cores = feature('numcores'); 
    % safe_workers = phys_cores - 2; % server: 14 local:8
    % 
    % if isempty(gcp('nocreate'))
    %     parpool('local', safe_workers); 
    % end
    for i = 1:num_frePoint
    % parfor i = 1:num_frePoint
        % slicing
        HB_i = HB_ctrl(:, :, i); RB = HB_i'*HB_i;
        HD_i = HD_ctrl(:, :, i);
        gammaD_i = 10;
        % eta_i = 0.5*sqrt(trace(RB)); % refer to [POTDC paper]
        eta_i = 1e-3*sqrt(trace(RB));
        wPRACC(:, i) = POTDC_solver(HB_i, HD_i, gammaD_i, eta_i, L, iterations);
    end
    % delete(gcp('nocreate')); % end the parallel pool
    w_opt.w = wPRACC;
    w_opt.scale = scaling_factor;
end

function w_opt = POTDC_solver(HB_i, HD_i, gammaD_i, eta_i, L, iterations)
    RB = HB_i' * HB_i;
    RD = HD_i' * HD_i; 
    [V, D] = eig(RB);
    [max_eigval_RB, idx] = max(diag(real(D)));    
    alpha_l = 1 / (1 - eta_i / sqrt(max(max_eigval_RB, eps)))^2;
    % alpha lower bound(加 eps 防止除零)
    max_eigval_RBD = real(max(eig((RD + gammaD_i * eye(L)) \ RB)));        
    k = 1.001 / (sqrt(max_eigval_RB) - eta_i);
    w0 = k * V(:, idx);  
    % construct feasible solution w0              
    alpha_u = max_eigval_RBD * real(w0' * (RD + gammaD_i * eye(L)) * w0);
    % upper bound
    alphaC = (alpha_u + alpha_l) / 2; % initial solutions
    W_current = zeros(L, L); 
    obj_prev = inf; % initial for object function    
    for j = 1:iterations  
        cvx_begin quiet 
            variable W(L, L) complex 
            variable alphaopt 
            
            minimize(real(trace((RD + gammaD_i * eye(L)) * W)))
            subject to
                real(trace(RB * W)) == alphaopt;
                W == semidefinite(L, L);
                alpha_l <= alphaopt;
                alphaopt <= alpha_u;
                eta_i^2 * real(trace(W)) + (sqrt(alphaC) - 1) + alphaopt * (1 / sqrt(alphaC) - 1) <= 0;
        cvx_end

        alphaC = alphaopt;
        W_current = W;        
        % convergence condition
        obj_curr = real(trace((RD + gammaD_i * eye(L)) * W));
        if j >= 2 && abs(obj_prev - obj_curr) <= 1e-12 % threshold
            break;
        end
        obj_prev = obj_curr;
    end   
    % --- recovery the rank==1 souliton ---
    [V_w, D_w] = eig(W_current);
    [lambda_max, max_idx_w] = max(diag(real(D_w)));
    v_max = V_w(:, max_idx_w);    
    w_opt = sqrt(max(lambda_max, 0)) * v_max;
end