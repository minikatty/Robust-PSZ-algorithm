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
    ATF_BZ_ctrl = operating_data.ATF_BZ.eval;

    target_indices = freq_params.target_indices;
    target_freqs = freq_params.target_freqs;
    num_target_freqs = length(target_indices);
    
    % --- 2. Pre-computation (if needed) ---
    % Example: Define target sound field for NSRE
    % This should be based on a virtual source at a specific location
    % For simplicity, let's assume it's the ATF of a single loudspeaker
    virtual_source_idx = eval_params.virtual_source_idx;
    p_d_all_freqs = ATF_BZ_ctrl(:, virtual_source_idx, :);
    
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