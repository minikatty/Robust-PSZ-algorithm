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
    eps_B_all = scaling_factor .* para.epsilon.B;
    eps_D_all = scaling_factor .* para.epsilon.D;
    rho = para.rho;
    mu = para.mu;
    alpha_AC = para.alpha;
    gamma = mu + rho * alpha_AC;    
    mosek_threads = [];
    rank_fraction_threshold = 0.98;
    normalization_multipliers = [0.1, 1, 5, 10, 0.01];
    capture_diagnostics = false;
    recovery_policy = 'adaptive'; % Keep historical production default.
    if isfield(para, 'racc_pm') && ...
            isfield(para.racc_pm, 'mosek_threads')
        mosek_threads = para.racc_pm.mosek_threads;
    end
    if isfield(para, 'racc_pm') && ...
            isfield(para.racc_pm, 'rank_fraction_threshold')
        rank_fraction_threshold = ...
            para.racc_pm.rank_fraction_threshold;
    end
    if isfield(para, 'racc_pm') && ...
            isfield(para.racc_pm, 'normalization_multipliers')
        normalization_multipliers = para.racc_pm.normalization_multipliers;
    end
    if isfield(para, 'racc_pm') && ...
            isfield(para.racc_pm, 'capture_diagnostics')
        capture_diagnostics = para.racc_pm.capture_diagnostics;
    end
    validateattributes(rank_fraction_threshold, {'numeric'}, ...
        {'scalar', 'real', 'finite', '>', 0, '<=', 1});
    validateattributes(normalization_multipliers, {'numeric'}, ...
        {'vector', 'real', 'finite', 'positive', 'nonempty'});
    validateattributes(capture_diagnostics, {'logical'}, {'scalar'});
    if isfield(para, 'racc_pm') && isfield(para.racc_pm, 'recovery_policy')
        recovery_policy = validatestring(para.racc_pm.recovery_policy, ...
            {'adaptive', 'principal'});
    end
    % Power constraint (inactive if large)
    e_w = 1e4; 
 
    wRACC_PM = zeros(L, num_frePoint);
    solution_status = strings(1, num_frePoint);
    selected_scale_multiplier = nan(1, num_frePoint);
    scale_attempt_status = cell(1, num_frePoint);
    % 1: default
    % 10, 100: enlarge
    % 0.1, 0.01: shrink
    % Defaults are unchanged; a singleton override permits controlled
    % normalization diagnostics without silently falling back to another.
    scale_multipliers = normalization_multipliers;
    solver_diagnostics = cell(1, num_frePoint);
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
    for i = loop_target(:).'
        eps_B = local_pick_value(eps_B_all, i);
        eps_D = local_pick_value(eps_D_all, i);
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
        attempt_statuses = strings(1, numel(scale_multipliers));
        w_tmp = zeros(L, 1);
        status_tmp = "Not run";
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
                                 'L', L, 'MB', MB, 'MD', MD, ...
                                 'mosek_threads', mosek_threads, ...
                                 'capture_diagnostics', capture_diagnostics, ...
                                 'recovery_policy', recovery_policy, ...
                                 'rank_fraction_threshold', ...
                                 rank_fraction_threshold);
        % Solving
            [w_tmp, status_tmp, diagnostics_tmp] = ...
                RACC_PM_LMI_solver(H_i, pB_d_i, epsilon_i, paraRACC_PM);
            if capture_diagnostics
                diagnostics_tmp.normalization_divisor = current_scale_factor;
                diagnostics_tmp.normalization_multiplier = scale_multipliers(retry);
                solver_diagnostics{i} = diagnostics_tmp;
            end
            attempt_statuses(retry) = string(status_tmp);

        % check status
            if strcmp(status_tmp, 'Solved') || strcmp(status_tmp, 'Inaccurate/Solved')
                % solved, save and break
                wRACC_PM(:, i) = w_tmp;
                solution_status(i) = string(status_tmp);
                selected_scale_multiplier(i) = scale_multipliers(retry);
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
        scale_attempt_status{i} = attempt_statuses;
    end    
    % delete(gcp('nocreate')); % end the parallel pool
    w_opt.w = wRACC_PM;
    w_opt.scale = scaling_factor;
    w_opt.status = solution_status;
    w_opt.normalization_multiplier = selected_scale_multiplier;
    w_opt.normalization_attempt_status = scale_attempt_status;
    w_opt.rank_fraction_threshold = rank_fraction_threshold;
    w_opt.recovery_policy = recovery_policy;
    if capture_diagnostics
        w_opt.solver_diagnostics = solver_diagnostics;
    end
end

function value = local_pick_value(source, frequency_index)
    if isscalar(source)
        value = source;
    elseif numel(source) >= frequency_index
        value = source(frequency_index);
    else
        error('RACC_PM_Sub:UncertaintyLengthMismatch', ...
            'No uncertainty radius is available for frequency index %d.', ...
            frequency_index);
    end
end


function [wopt_RACC_PM, status, diagnostics] = ...
        RACC_PM_LMI_solver(H_i, pB_d, epsilon_i, paraRACC_PM)
    % Unpack parameters
    rho = paraRACC_PM.rho; 
    gamma_i = paraRACC_PM.gamma_i; 
    ew = paraRACC_PM.ew; 
    L = paraRACC_PM.L; 
    MB = paraRACC_PM.MB; 
    MD = paraRACC_PM.MD;
    mosek_threads = paraRACC_PM.mosek_threads;
    rank_fraction_threshold = ...
        paraRACC_PM.rank_fraction_threshold;
  
    HB_i = H_i.B;
    HD_i = H_i.D;
    epsB = epsilon_i.B;
    epsD = epsilon_i.D;
    
    % Noise floor for stability (optional constant term)
    sigma_sq = 1e-8; 

    cvx_begin sdp quiet
        % Recommended solver for decomposed problems
        cvx_solver mosek 
        if ~isempty(mosek_threads)
            cvx_solver_settings('MSK_IPAR_NUM_THREADS', mosek_threads);
        end
        
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
                % vec(DeltaH) uses a non-conjugating row-block transpose.
                u_vec_B = ((1 - rho) * h_m * W - p_m * w').';
                
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
                u_vec_D = (gamma_i * h_k * W).';
                
                % The LMI Block
                [ tau_D * eye(L) - gamma_i * W.',   u_vec_D;
                  u_vec_D',                         delta_D(k) ] >= 0;
            end
            
    cvx_end
    
    status = cvx_status;    
    diagnostics = struct();
    if paraRACC_PM.capture_diagnostics
        diagnostics = struct('status', status, 'objective_scaled', cvx_optval, ...
            'sigma_sq_scaled', sigma_sq, 'W_tilde', W_tilde, ...
            'tB', tB, 'tD', tD, 'tau_B', tau_B, 'tau_D', tau_D, ...
            'delta_B', delta_B, 'delta_D', delta_D);
    end
    % --- Post-Processing / Rank-1 Approximation ---
    if ~strcmp(status, 'Solved') && ~strcmp(status, 'Inaccurate/Solved')
        wopt_RACC_PM = zeros(L, 1);
        % warning('CVX problem not solved: %s', status);
        return;
    end
    
    % Check whether the dominant eigencomponent captures the lifted matrix.
    % This full-spectrum criterion avoids classifying a nearly rank-one
    % solution from lambda_1/lambda_2 alone.
    W_tilde = (W_tilde + W_tilde') / 2;
    [V, eig_vals] = eig(full(W_tilde), 'vector');
    eig_vals = max(real(eig_vals), 0);
    [eig_vals, eig_order] = sort(eig_vals, 'descend');
    V = V(:, eig_order);
    lambda1 = eig_vals(1);
    rank_fraction = lambda1 / max(sum(eig_vals), realmin);
    lifted_mean = W_tilde(1:L, L+1);

    if abs(V(L+1, 1)) > 1e-10
        w_extracted = V(1:L, 1) / V(L+1, 1);
    else
        % Stable fallback if the homogeneous coordinate is degenerate.
        w_extracted = lifted_mean;
    end
    if norm(w_extracted)^2 > ew
        w_extracted = w_extracted * sqrt(ew) / norm(w_extracted);
    end
    
    if strcmp(paraRACC_PM.recovery_policy, 'principal') || ...
            rank_fraction >= rank_fraction_threshold
        % Effectively rank one: deterministic principal recovery.
        wopt_RACC_PM = w_extracted;
        recovery_method = 'Principal eigenvector';
    else
        % Protective Gaussian randomization. For the affine lifting
        % [W,w;w^H,1], sample candidates that preserve both moments:
        % E[w_k]=w and E[w_k*w_k^H]=W.
        num_rand = 1000;
        lifted_covariance = W_tilde(1:L, 1:L) - ...
            lifted_mean * lifted_mean';
        lifted_covariance = ...
            (lifted_covariance + lifted_covariance') / 2;
        [U, covariance_eigenvalues] = eig( ...
            full(lifted_covariance), 'vector');
        covariance_eigenvalues = max(real(covariance_eigenvalues), 0);
        covariance_square_root = U * ...
            diag(sqrt(covariance_eigenvalues));

        % Always retain the deterministic candidate as a fallback, so the
        % randomization branch cannot return a worse nominal proxy.
        Hw_B = HB_i * w_extracted;
        err_B = pB_d - Hw_B;
        Hw_D = HD_i * w_extracted;
        best_obj = real(err_B' * err_B - rho * (Hw_B' * Hw_B) + ...
            gamma_i * (Hw_D' * Hw_D));
        w_best = w_extracted;
        
        for k = 1:num_rand
            r = (randn(L, 1) + 1j*randn(L, 1)) / sqrt(2);
            w_k = lifted_mean + covariance_square_root * r;
            
            % Enforce Power
            pwr = norm(w_k)^2;
            if pwr > ew
                w_k = w_k * sqrt(ew / pwr);
            end
            
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
        recovery_method = 'Gaussian randomization';
    end
    if paraRACC_PM.capture_diagnostics
        diagnostics.eigenvalues = eig_vals;
        diagnostics.principal_anchor_valid = abs(V(L+1, 1)) > 1e-10;
        diagnostics.eigenvalue_ratio = lambda1 / max(eig_vals(2), realmin);
        diagnostics.relative_rank_one_error = norm(W_tilde - ...
            lambda1*V(:, 1)*V(:, 1)', 'fro') / ...
            max(norm(W_tilde, 'fro'), realmin);
        diagnostics.rank_fraction = rank_fraction;
        diagnostics.recovery_method = recovery_method;
        diagnostics.principal_filter = w_extracted;
        diagnostics.lifted_mean = lifted_mean;
        diagnostics.recovered_filter = wopt_RACC_PM;
    end
end
