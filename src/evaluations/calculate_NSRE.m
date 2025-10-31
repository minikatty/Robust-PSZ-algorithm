function NSRE = calculate_NSRE(w, H_B, p_d)
% CALCULATE_NSRE - Computes the Normalized Sound Reproduction Error (NSRE) in dB.
% the desired virtual source field: p_d
    p_reproduced = H_B * w;
    
    % Add a small epsilon to prevent division by zero if target field is silent
    epsilon = 1e-12;
    
    error_power = norm(p_reproduced - p_d)^2;
    target_power = norm(p_d)^2;
    
    NSRE = 10 * log10(error_power / (target_power + epsilon) + epsilon);
end