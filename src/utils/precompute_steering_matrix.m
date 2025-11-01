function S = precompute_steering_matrix(mic_positions, freqs, c)
% PRECOMPUTE_PLANARITY_MATRIX - Precomputes the S matrix for planarity calculation.
    BZ_origin_shift = [2.6,2,0];
    planarity_eval_mic_position(array);
    num_mics = size(mic_positions, 1);
    num_freqs = length(freqs);
    mic_positions = mic_positions - BZ_origin_shift;
    
    S = zeros(360, num_mics, num_freqs);
    angles_rad = (0:359)' * (2*pi/360);
    ui_all = [sin(angles_rad), cos(angles_rad)];
    
    for f_idx = 1:num_freqs
        k = 2 * pi * freqs(f_idx) / c; % freq2wavenumber
        % Vectorized calculation for a single frequency
        S(:, :, f_idx) = exp(-1i * k * (ui_all * mic_positions(:,1:2)')) / num_mics;
    end
end