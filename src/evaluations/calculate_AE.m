function AE = calculate_AE(w, H_B, virtual_source_idx)
% CALCULATE_AE - Computes the Array Effort (AE) in dB.
    
    w0 = zeros(size(w));
    w0(virtual_source_idx) = 1;
    
    % Add a small epsilon for numerical stability
    epsilon = 1e-12;
    
    effort_ratio = (w' * w);
    RB = (H_B' * H_B);
    energy_ratio = (w0' * RB * w0) / (w' * RB * w + epsilon);
    
    AE = 10 * log10(effort_ratio * energy_ratio + epsilon);
end