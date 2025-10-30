function performance = evaluate_performance(filters_w, operating_data, freq_params)
% EVALUATE_PERFORMANCE - Calculates various performance metrics for a given set of filters.

    % Extract the necessary ATFs for evaluation (usually evaluation points)
    ATF_BZ_eval = operating_data.ATF_BZ.eval;
    ATF_DZ_eval = operating_data.ATF_DZ.eval;
    
    % Extract target frequency indices
    target_indices = freq_params.target_indices;
    
    % Select the ATFs at the target frequencies
    H_B = ATF_BZ_eval(:, :, target_indices); % Bright Zone ATF
    H_D = ATF_DZ_eval(:, :, target_indices); % Dark Zone ATF
    
    % Reshape filters_w if it's for multiple frequencies
    % w = filters_w(:, target_indices); % This depends on your filter format
    w = filters_w; % Assuming filters_w is already [N_src x N_freq]
    
    % --- Calculate Sound Field at each target frequency ---
    % This is a simplified example, your calculation might be more complex
    % p_B = H_B * w; % This needs to be a proper matrix multiplication for each freq
    % p_D = H_D * w;
    
    % Placeholder calculations - YOU MUST REPLACE THESE
    avg_energy_B = 10; % Replace with actual calculation, e.g., mean(abs(p_B).^2)
    avg_energy_D = 0.1;

    % --- Calculate Metrics ---
    % Metric 1: Acoustic Contrast
    performance.AcousticContrast = 10 * log10(avg_energy_B / avg_energy_D);
    
    % Metric 2: Quiet Zone SPL (assuming Dark Zone is the quiet zone)
    % Assuming a reference pressure of 20e-6 Pa for SPL
    p_ref = 20e-6;
    performance.QuietZoneSPL = 10 * log10(avg_energy_D / p_ref^2);
    
    % Metric 3: Control Effort
    performance.ControlEffort = sum(abs(w(:)).^2);
    
end


% src/evaluation/evaluate_performance.m

function performance = evaluate_performance(filters_w, operating_data, freq_params, eval_params)
% EVALUATE_PERFORMANCE - Calculates all performance metrics for a given set of filters.
%
% Inputs:
%   filters_w     - Control filters [N_src x N_freq_target]
%   operating_data- Struct containing ATFs for the operating condition
%   freq_params   - Struct with frequency parameters
%   eval_params   - Struct with evaluation-specific parameters (e.g., virtual source)

    % --- 1. Extract Data ---
    ATF_BZ_eval = operating_data.ATF_BZ.eval;
    ATF_DZ_eval = operating_data.ATF_DZ.eval;
    
    target_indices = freq_params.target_indices;
    target_freqs = freq_params.target_freqs;
    num_target_freqs = length(target_indices);
    
    % --- 2. Pre-computation (if needed) ---
    % Example: Define target sound field for NSRE
    % This should be based on a virtual source at a specific location
    % For simplicity, let's assume it's the ATF of a single loudspeaker
    virtual_source_idx = eval_params.virtual_source_idx;
    p_d_all_freqs = ATF_BZ_eval(:, virtual_source_idx, :);
    
    % Example: Precompute S matrix for Planarity
    % This is computationally expensive and should ideally be done once outside the loop
    if ~isfield(eval_params, 'S_matrix')
        error('S_matrix for Planarity is not provided in eval_params.');
    end
    S_matrix_all_freqs = eval_params.S_matrix;
    
    % --- 3. Initialize result vectors ---
    AC = zeros(num_target_freqs, 1);
    NSRE = zeros(num_target_freqs, 1);
    AE = zeros(num_target_freqs, 1);
    Planarity = zeros(num_target_freqs, 1);

    % --- 4. Loop through each target frequency and calculate metrics ---
    for i = 1:num_target_freqs
        f_idx = target_indices(i);
        
        % Extract data for the current frequency
        w_f = filters_w(:, i);
        H_B_f = ATF_BZ_eval(:, :, f_idx);
        H_D_f = ATF_DZ_eval(:, :, f_idx);
        p_d_f = p_d_all_freqs(:, 1, f_idx);
        S_f = S_matrix_all_freqs(:, :, f_idx);
        
        % Calculate metrics using the helper functions
        AC(i) = calculate_ac(w_f, H_B_f, H_D_f);
        NSRE(i) = calculate_nsre(w_f, H_B_f, p_d_f);
        AE(i) = calculate_ae(w_f, H_B_f, virtual_source_idx);
        Planarity(i) = calculate_planarity(w_f, H_B_f, S_f);
    end
    
    % --- 5. Package into output structure ---
    % We can return the full curves or the average over the frequency band
    performance.AC = mean(AC);
    performance.NSRE = mean(NSRE);
    performance.AE = mean(AE);
    performance.Planarity = mean(Planarity);
    
    % Optional: return the full curves as well
    performance.AC_curve = AC;
    performance.freq_axis_target = target_freqs;
end