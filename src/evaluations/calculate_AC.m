function AC = calculate_AC(w, H_B, H_D)
% CALCULATE_AC - Computes the Acoustic Contrast (AC) in dB.
    M_B = size(H_B, 1);
    M_D = size(H_D, 1);
    
    numerator = (w' * (H_B' * H_B) * w);
    denominator = (w' * (H_D' * H_D) * w);
    
    % Add a small epsilon to prevent log10(0) or division by zero
    epsilon = 1e-8;
    
    mu = (M_D * numerator) / (M_B * denominator + epsilon);
    AC = 10 * log10(mu + epsilon); % dB
end