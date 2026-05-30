function [w_opt] = RACC_PM_Sub(HB_ctrl, HD_ctrl, H_desired, para, fre_idces)
% Decompose the large scale problem into sub-contraints
% Optimization Problem Formulation (Eq. 50a - 50h)
% Uses Schur Complement Decomposition for efficient solving.
% ------------------------------------------------
% Variables:
%     W_tilde: (L+1)x(L+1) Hermitian matrix
%     t_B, t_D: Scalar auxiliary variables (objectives)
%     tau_B, tau_D: Scalar S-procedure variables
%     delta_B, delta_D: Slack variable vectors (size M_B and M_D)
% 
% Objective:
%     (50a) min (t_B + t_D)
% 
% Constraints:
%     -- Bright Zone (BZ) --
%     (50b) Robust LMIs: For m = 1..M_B:
%           [ tau_B*I - (1-rho)*W.T   u_Bm  ]
%           [ u_Bm'                   delta ]  >= 0
%     (50c) Slack Sum: sum(delta_B) <= -c_B - tau_B * eps_B^2
% 
%     -- Dark Zone (DZ) --
%     (50d) Robust LMIs: For k = 1..M_D:
%           [ tau_D*I - gamma*W.T     u_Dk  ]
%           [ u_Dk'                   delta ]  >= 0
%     (50e) Slack Sum: sum(delta_D) <= -c_D - tau_D * eps_D^2
% 
%     -- Structural & Power --
%     (50f) PSD: W_tilde >= 0, W_tilde[L+1, L+1] == 1
%     (50g) Power: Tr(Sigma * W_tilde) <= e_w
%     (50h) Non-negative: tau_B >= 0, tau_D >= 0
%
% Inputs:
%   HB_ctrl : (MB x L x Freq)
%   HD_ctrl : (MD x L x Freq)
%   H_desired: (MB x Freq)
%   para    : Structure containing epsilon, rho, mu, alpha, etc.
%
% Outputs:
%   w_opt   : (L x 1) Optimal filter coefficients
%   status  : CVX status string

    [MB, L, num_frePoint] = size(HB_ctrl);
    [MD, ~] = size(HD_ctrl);

    if nargin < 5 || isempty(fre_idces)
        loop_target = 1:num_frePoint; 
    else
        loop_target = fre_idces;          
    end
    
    scaling_factor = para.scale; % control the bound of uncertainty
    % eps_B = 0.5 * sqrt(trace(RB)); % scaling_factor .* para.epsilon.B
    % eps_D = 0.5 * sqrt(trace(RD)); % scaling_factor .* para.epsilon.D
    rho = para.rho;
    mu = para.mu;
    alpha_AC = para.alpha;
    gamma = mu + rho * alpha_AC;    
    % Power constraint (inactive if large)
    e_w = 1e4; 
 
    wRACC_PM = zeros(L, num_frePoint);
    solution_status =  strings(1, num_frePoint);
    % 1: default
    % 10, 100: enlarge
    % 0.1, 0.01: shrink
    scale_multipliers =[0.1, 1, 5, 10, 0.01]; 
    cvx_clear;

    % Optional: Parallel Pool Check but fail in large scale problem
    % if isempty(gcp('nocreate')) && num_frePoint > 1
    %     parpool('local'); 
    % end
    % fprintf('Starting Robust Optimization (Decomposed Form)...\n');

    % phys_cores = feature('numcores'); 
    % if phys_cores >= 12
    %     safe_workers = 12;   % server: 14
    % else
    %     safe_workers = 4;   % local:8
    % end
    % 
    % if isempty(gcp('nocreate'))
    %     parpool('local', safe_workers); 
    % end
    

    % parfor i = 1: num_frePoint
    for i = loop_target
        HB = squeeze(HB_ctrl(:, :, i));
        HD = squeeze(HD_ctrl(:, :, i));
        RB = HB'*HB;
        RD = HD'*HD;  
        eps_B = 1e-2 * sqrt(norm(RB,'fro')); % scaling_factor .* para.epsilon.B:5e-3
        eps_D = 1e-2 * sqrt(norm(RD,'fro')); % scaling_factor .* para.epsilon.D 
        % 因为wcACC曾说不考虑暗区的不确定性反而有助于提高AC指标，所以这里尝试暗区用很小的不确定性
        %0.001 * sqrt(trace(RD))
        HB_raw = HB_ctrl(:, :, i);
        HD_raw = HD_ctrl(:, :, i);
        pB_raw = H_desired(:, i); % raw pressure
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


function [wopt_RACC_PM, status] = RACC_PM_LMI_solver(H_i, pB_d, epsilon_i, paraRACC_PM)
    % Unpack parameters
    rho = paraRACC_PM.rho; 
    gamma_i = paraRACC_PM.gamma_i; 
    ew = paraRACC_PM.ew; 
    L = paraRACC_PM.L; 
    MB = paraRACC_PM.MB; 
    MD = paraRACC_PM.MD;
  
    HB_i = H_i.B;
    HD_i = H_i.D;
    epsB = epsilon_i.B;
    epsD = epsilon_i.D;
    
    % Noise floor for stability (optional constant term)
    sigma_sq = 1e-8; 

    cvx_begin sdp quiet
        % Recommended solver for decomposed problems
        cvx_solver mosek 
        
        % --- Variables ---
        % Augmented matrix W_tilde = [W, w; w', 1]
        variable W_tilde(L+1, L+1) hermitian
        
        % Auxiliary variables for objective
        variable tB
        variable tD
        
        % S-Procedure multipliers (scalars)
        variable tau_B nonnegative
        variable tau_D nonnegative
        
        % Decomposed Slack Variables 
        variable delta_B(MB)
        variable delta_D(MD)
        
        % --- Objective ---
        minimize( tB + tD )
        
        subject to
            % --- Weight Matrix Structural Constraints ---
            W_tilde >= 0;
            W_tilde(L+1, L+1) == 1;
            
            % Extract W and w using expressions 
            W = W_tilde(1:L, 1:L);
            w = W_tilde(1:L, L+1);
            
            % Power Constraint
            real(trace(W)) <= ew;
            
            % =========================================================
            % 1. Robust Constraint for Bright Zone (Decomposed)
            % =========================================================
            
            % Calculate Constant c_B
            % c_B = (1-rho)Tr(H W H') - 2Re(p' H w) + p'p + sigma^2 - tB
            % Note: trace(H*W*H') can be computed efficiently via sum of quadratic forms
            term_quad_B = real(trace(HB_i * W * HB_i')); 
            term_lin_B  = 2 * real(pB_d' * HB_i * w);
            const_term_p = pB_d' * pB_d;
            cB = (1 - rho) * term_quad_B - term_lin_B + const_term_p + sigma_sq - tB;
            % (a) Summation Constraint
            sum(delta_B) <= -cB - tau_B * (epsB^2);
            
            % (b) Parallel Block Constraints
            for m = 1:MB
                % Data for m-th constraint
                h_m = HB_i(m, :);   % 1xL vector
                p_m = pB_d(m);      % Scalar
                
                % Construct vector u_{B,m}
                % Formula: u_{B,m} = [ (1-rho)*h_m*W - p_m*w^H ]^H
                % Note: w' in MATLAB is conjugate transpose (w^H).
                % We need w^T here? No, p_m*w' is scalar product in derivation.
                % Let's follow the vector form derived:
                % u_vec = ( (1-rho)*h_m*W - p_m*w' )' -> Transforms row to col
                u_vec_B = ((1 - rho) * h_m * W - p_m * w')';
                
                % The LMI Block:
                % Matrix Omega_B = tau_B*I - (1-rho)*W^T
                % Note: Since W is Hermitian, W^T = conj(W). 
                % In CVX, W.' is transpose (non-conjugate). 
                % If W is real, W.'=W. If complex, W.' != W'.
                % Correct mapping: The math derivation W^T corresponds to conj(W) if W Hermitian.
                % Using W.' (Transpose) creates the correct structure for [conj(W) u; u' delta].
                
                [ tau_B * eye(L) - (1-rho) * W.',   u_vec_B;
                  u_vec_B',                         delta_B(m) ] >= 0;
            end
            
            % =========================================================
            % 2. Robust Constraint for Dark Zone (Decomposed)
            % =========================================================
            
            % Calculate Constant c_D
            % c_D = gamma * Tr(H W H') - tD
            cD = gamma_i * real(trace(HD_i * W * HD_i')) - tD;
            
            % (a) Summation Constraint
            sum(delta_D) <= -cD - tau_D * (epsD^2);
            
            % (b) Parallel Block Constraints
            for k = 1:MD
                h_k = HD_i(k, :);
                
                % Construct vector u_{D,k}
                % Formula: u_{D,k} = [ gamma * h_k * W ]^H
                u_vec_D = (gamma_i * h_k * W)';
                
                % The LMI Block
                [ tau_D * eye(L) - gamma_i * W.',   u_vec_D;
                  u_vec_D',                         delta_D(k) ] >= 0;
            end
            
    cvx_end
    
    status = cvx_status;    
    % --- Post-Processing / Rank-1 Approximation ---
    if ~strcmp(status, 'Solved') && ~strcmp(status, 'Inaccurate/Solved')
        wopt_RACC_PM = zeros(L, 1);
        % warning('CVX problem not solved: %s', status);
        return;
    end
    
    % Check Rank
    eig_vals = sort(eig(W_tilde), 'descend');
    % Handle numerical zeros
    lambda1 = abs(eig_vals(1));
    lambda2 = abs(eig_vals(2));
    if lambda2 < 1e-10, lambda2 = 1e-10; end
    rank_ratio = lambda1 / lambda2;
    
    if rank_ratio > 1e3 % Considered Rank-1
        % Extract principal eigenvector
        [V, ~] = eigs(W_tilde, 1);
        % Normalize to satisfy W_tilde(L+1,L+1)=1
        w_extracted = V(1:L) / V(L+1);
        
        % Enforce power constraint
        if norm(w_extracted)^2 > ew
             w_extracted = w_extracted * sqrt(ew) / norm(w_extracted);
        end
        wopt_RACC_PM = w_extracted;
    else
        % Gaussian Randomization
        % fprintf('High rank (ratio %.1f). Randomizing...\n', rank_ratio);
        num_rand = 1000;
        [U, S] = eig(W_tilde);
        best_obj = inf;
        w_best = zeros(L, 1);
        
        % Pre-compute metric constants to speed up loop
        % We use the WORST-CASE metric proxy (Nominal + Penalty)
        
        for k = 1:num_rand
            % Generate candidate w
            r = (randn(L+1, 1) + 1j*randn(L+1, 1)) / sqrt(2);
            w_cand_aug = U * sqrt(S) * r;
            w_k = w_cand_aug(1:L);
            
            % Enforce Power
            pwr = norm(w_k)^2;
            if pwr > ew
                w_k = w_k * sqrt(ew / pwr);
            end
            
            % Evaluate Objective (Proxy: Nominal Cost)
            % f_B = ||p - Hw||^2 - rho*||Hw||^2
            % f_D = gamma * ||Hw||^2
            
            Hw_B = HB_i * w_k;
            err_B = pB_d - Hw_B;
            fB_val = (err_B' * err_B) - rho * (Hw_B' * Hw_B);
            
            Hw_D = HD_i * w_k;
            fD_val = gamma_i * (Hw_D' * Hw_D);
            
            current_obj = real(fB_val + fD_val);
            
            if current_obj < best_obj
                best_obj = current_obj;
                w_best = w_k;
            end
        end
        wopt_RACC_PM = w_best;
    end
end