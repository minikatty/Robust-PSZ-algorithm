function [w_opt] = RACC_PM_GLS(HB_ctrl, HD_ctrl, H_desired, para, fre_idces)

%ROBUST_SFC_DESIGN Solves the robust sound field control problem via SDR.
%  global large scale (GLS)
%   Minimize  t_B + t_D
%   s.t.      Robust Bright Zone Constraint (via S-Procedure LMI)
%             Robust Dark Zone Constraint (via S-Procedure LMI)
%             AE Constraint
%
% Inputs:
%   H_B    : (Mb x L x Freq) Nominal ATF matrix for Bright Zone
%   H_D    : (Md x L x Freq) Nominal ATF matrix for Dark Zone
%   p_d    : (Mb x 1 x Freq) Desired sound pressure vector
%   par    : Structure containing parameters:
%            para.eps_B  - Uncertainty bound for H_B
%            para.eps_D  - Uncertainty bound for H_D
%            par.rho    - Weighting factor
%            para.mu     - Parameter for Dark Zone
%            para.alpha  - Parameter for Dark Zone
%            para.num_rand - Number of Gaussian randomizations (e.g., 1000)
%
% Outputs:
%   w_opt   : (L x 1) Optimal filter coefficients
%   status  : CVX status string

    [MB, L, num_frePoint] = size(HB_ctrl);
    [MD, ~, ~] = size(HD_ctrl);

    if nargin < 5 || isempty(fre_idces)
        loop_target = 1:num_frePoint; 
    else
        loop_target = fre_idces;          
    end

    scaling_factor = para.scale; % control the bound of uncertainty
    % eps_B = scaling_factor.* para.epsilon.B;
    % eps_D = scaling_factor.* para.epsilon.D;
    rho = para.rho;
    mu = para.mu;
    alpha_AC = para.alpha;
    gamma = mu + rho * alpha_AC;
    % ew    = para.ew; % AE power constraint
    e_w = 1e4; % a large value keeps the contraint inactive 
   
    wRACC_PM = zeros(L, num_frePoint, 'like', 1i);
    solution_status = strings(1, num_frePoint);
    % 1: default
    % 10, 100: enlarge
    % 0.1, 0.01: shrink
    scale_multipliers =[0.1, 1, 5, 10, 0.01]; 
    cvx_clear;

    % phys_cores = feature('numcores'); 
    % 
    % if phys_cores >= 12
    %     safe_workers = 6;   % server: 14
    % else
    %     safe_workers = 1;   % local:8
    % end
    % 
    % if isempty(gcp('nocreate'))
    %     parpool('local', safe_workers); 
    % end

    % parfor i = 1:num_frePoint
    for i = loop_target
        HB = squeeze(HB_ctrl(:, :, i));
        HD = squeeze(HD_ctrl(:, :, i));
        RB = HB'*HB;
        RD = HD'*HD;
        eps_B = 1e-2 * sqrt(norm(RB,'fro')); % scaling_factor .* para.epsilon.B
        eps_D = 1e-2 * sqrt(norm(RD,'fro')); % scaling_factor .* para.epsilon.D
        HB_raw = HB_ctrl(:, :, i);
        HD_raw = HD_ctrl(:, :, i);
        pB_raw = H_desired(:, i); % raw pressure
        % to make rho * ||H||^2  ≈ 1
        H_energy = norm([HB_raw; HD_raw], 'fro');            
        if H_energy < 1e-10, H_energy = 1; end 
        % aviod anormal small number        
        % alpha_i = ||H|| * sqrt(rho)
        base_scale_factor = H_energy * sqrt(rho);
        % Prepare data for current frequency
        % Scaling        
        for retry = 1:length(scale_multipliers)
            current_scale_factor = base_scale_factor * scale_multipliers(retry);
            % Prepare data for current frequency            
            H_i = struct();
            H_i.B = HB_raw / current_scale_factor;
            H_i.D = HD_raw / current_scale_factor;
            pB_d_i = pB_raw / current_scale_factor;
    
            epsilon_i = struct('B', eps_B / current_scale_factor, ...
                               'D', eps_D / current_scale_factor);   
            gamma_i = gamma(i);
            
            paraRACC_PM = struct('rho', rho, 'gamma_i', gamma_i, 'ew', e_w, ...
                                 'L', L, 'MB', MB, 'MD', MD);        
        % Solving
            [w_tmp, status_tmp] = RACC_PM_LMI_solver(H_i, pB_d_i, epsilon_i, paraRACC_PM);

        % check status
            if strcmp(status_tmp, 'Solved') || strcmp(status_tmp, 'Inaccurate/Solved')
                % solved, save and break
                wRACC_PM(:, i) = w_tmp;
                solution_status(i) = status_tmp;                
                % if try out works, print message
                if retry > 1
                    fprintf('-> Freqency idx %d-th: successful after adjusting scaling facor: %g \n', i, scale_multipliers(retry));
                end
                break; 
            end
            % all tryout failed, print the frequency bin
            if retry == length(scale_multipliers)
                wRACC_PM(:, i) = w_tmp; 
                solution_status(i) = status_tmp;
                warning('frequency %d: try %d factors，all failed! Terminal status: %s', i, retry, status_tmp);
            end        
        % Simple progress print
        if mod(i, 10) == 0 || i == 1
            fprintf('Processed Frequency %d / %d\n', i, num_frePoint);
        end
        end
    end    
    % delete(gcp('nocreate')); % end the parallel pool
    w_opt.w = wRACC_PM;
    w_opt.scale = scaling_factor;
    w_opt.status = solution_status;
end

function [wopt_RACC_PM,status] = RACC_PM_LMI_solver(H_i, pB_d, epsilon_i, paraRACC_PM)
    % HB_i (亮区传递函数), HD_i (暗区传递函数), pB_d (明区目标声场)
    % mu (权衡因子), epsilon (扰动量参数), ew (能量约束上限)
    % L (扬声器数量)        
    rho = paraRACC_PM.rho; gamma_i = paraRACC_PM.gamma_i; ew = paraRACC_PM.ew; 
    L = paraRACC_PM.L; MB = paraRACC_PM.MB; MD = paraRACC_PM.MD;
    I_MB = eye(MB); % MB×MB
    I_MD = eye(MD); % Identity matrices for Kronecker products
    HB_i = H_i.B;
    HD_i = H_i.D;
    epsilonB_i = epsilon_i.B;
    epsilonD_i = epsilon_i.D;    

    cvx_begin sdp quiet
        cvx_solver sdpt3
        % cvx_solver mosek 
        % --- Variables ---
        % Augmented matrix W_tilde = [W, w; w', 1]
        variable W_tilde(L+1, L+1) hermitian
        variable tB
        variable tD
        variable lambda_B nonnegative
        variable lambda_D nonnegative            
        % Extract sub-blocks for cleaner expressions
        % Note: Using expressions to avoid new variables
        expression W(L, L)
        expression w(L, 1)
        expression Q_B(MB*L, MB*L)
        expression Q_D(MD*L, MD*L)
        W = W_tilde(1:L, 1:L);
        w = W_tilde(1:L, L+1);    
        % Objective function
        minimize( tB + tD )    
        subject to
            % --- Structural Constraints ---
            W_tilde >= 0;
            W_tilde(L+1, L+1) == 1;                
            % AE Constraint: Tr(W) <= ew % not consider for now <==> a large ew
            real(trace(W)) <= ew; 
            % --- Robust Constraint for BZ -- 
            % Q_B = (1-rho) * (W^T kron I_Mb)               
            Q_B = (1 - rho) * kron(conj(W), I_MB); 
            % cvx's bug
            % Q_B size: (MB*L) × (MB*L)            
            % for i = 1:L
            %     rows = (i-1)*MB + (1:MB);
            %     for j = 1:L
            %         cols = (j-1)*MB + (1:MB);
            %         Q_B(rows, cols) = (1-rho) * W(j,i) * I_MB;   %  W(j,i)=W^T
            %     end
            % end

            % u_B = vec( (1-rho)H_B*W - p_d*w' )
            u_B = (1 - rho) * vec(HB_i * W - pB_d * w');              

            % c_B = (1-rho)Tr(H W H') - 2Re(p' H w) + p'p - tB
            cB = (1 - rho) * real(trace(HB_i * W * HB_i')) ...
                - 2 * real(pB_d' * HB_i * w) + (pB_d' * pB_d) - tB;

            % LMI for f_B
            [ lambda_B * eye(MB*L) - Q_B,   -u_B;
              -u_B',                        -cB - lambda_B * epsilonB_i^2 ] >= 0;

            % --- Robust Constraint for DZ ---
            % Q_D = gamma_D * (W^T kron I_Md)
            Q_D = gamma_i * kron(conj(W), I_MD);

            % u_D = vec( gamma_D * H_D * W )
            u_D = vec(gamma_i * HD_i * W);                
            % c_D calculation
            cD = gamma_i * real(trace(HD_i * W * HD_i')) - tD;                
            % LMI for Eve
            [ lambda_D * eye(MD*L) - Q_D,   -u_D;
              -u_D',                        -cD - lambda_D * epsilonD_i^2 ] >= 0;              
    cvx_end
    status = cvx_status;
    % opt_val = cvx_optval; % the optval of objective
    if ~strcmp(status, 'Solved') && ~strcmp(status, 'Inaccurate/Solved')
        wopt_RACC_PM = zeros(L, 1, 'like', 1i);
        warning('CVX problem not solved: %s', status);
        return;
    end
    % Rank-1 Recovery (Gaussian Randomization)
    eig_vals = sort(eig(W_tilde), 'descend');
    rank_ratio = abs(eig_vals(1)) / (abs(eig_vals(2)) + 1e-10);
    if rank_ratio > 1e3
    % --- Rank-1 (Tight) ---
    % Extract principal eigenvector directly, Normalize by the last element to ensure [w; 1] structure
        [V, ~] = eigs(W_tilde, 1);
        wopt_RACC_PM = V(1:L) / V(L+1);        
    % Ensure power constraint (numerical safety)
        if norm(wopt_RACC_PM)^2 > ew
            wopt_RACC_PM = wopt_RACC_PM * sqrt(ew) / norm(wopt_RACC_PM);
        end        
    else
    % --- High Rank (Gaussian Randomization) ---
    % fprintf('High rank solution (ratio %.1f). Running randomization...\n', rank_ratio);        
        num_rand = 2000;
        [U, S] = eig(W_tilde);            
        best_obj = inf;
        w_best = zeros(L, 1, 'like', 1i);
       
        for k = 1:num_rand
            % Generate candidate
            % random direction
            r = (randn(L+1, 1) + 1j*randn(L+1, 1)) / sqrt(2);
            w_cand_aug = U * sqrt(S) * r;            
            % Extract w
            w_k = w_cand_aug(1:L); % Simple truncation usually works best for SDR            
            % Scale to satisfy Power Constraint
            AE_pwr = norm(w_k)^2;
            if AE_pwr > ew
                w_k = w_k * sqrt(ew / AE_pwr);
            end

            % Evaluate Metric
            % Strictly speaking, we should evaluate the WORST-CASE cost here.
            % However, solving a maximization problem for each candidate is slow.
            % We use the Nominal Cost + Penalty as a proxy metric for selection.

            % Nominal Bright Zone Cost
            err_B = pB_d - HB_i * w_k;
            fB_nom = norm(err_B)^2 - rho * norm(HB_i * w_k)^2;

            % Nominal Dark Zone Cost
            fD_nom = gamma_i * norm(HD_i * w_k)^2;

            % Proxy Objective
            current_obj = fB_nom + fD_nom;

            if current_obj < best_obj
                best_obj = current_obj;
                w_best = w_k;
            end
        end
            wopt_RACC_PM = w_best;
    end
    % %% for debug
    % status = 'Solved';
    % wopt_RACC_PM = zeros(48,1,'like',1+1i);
end