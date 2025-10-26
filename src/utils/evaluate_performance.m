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