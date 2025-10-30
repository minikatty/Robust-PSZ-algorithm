function S = precompute_planarity_matrix(mic_positions, freqs, c)
% PRECOMPUTE_PLANARITY_MATRIX - Precomputes the S matrix for planarity calculation.
    
    num_mics = size(mic_positions, 1);
    num_freqs = length(freqs);
    
    S = zeros(360, num_mics, num_freqs);
    angles_rad = (0:359)' * (2*pi/360);
    ui_all = [sin(angles_rad), cos(angles_rad)];
    
    for f_idx = 1:num_freqs
        k = 2 * pi * freqs(f_idx) / c;
        % Vectorized calculation for a single frequency
        S(:, :, f_idx) = exp(-1i * k * (ui_all * mic_positions(:,1:2)')) / num_mics;
    end
end