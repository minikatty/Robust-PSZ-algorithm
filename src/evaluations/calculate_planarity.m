function planarity = calculate_planarity(w, H_B, S_matrix)
% CALCULATE_PLANARITY - Computes the planarity of the sound field.
% H_B = [M,L,freq]
% w = [L,freq]
    p_B = H_B * w;
    psi = S_matrix * p_B;
    
    flux = 0.5 * abs(psi).^2; % Energy at each angle
    
    [~, max_idx] = max(flux);
    
    % Angles vector (0 to 359 degrees)
    angles_rad = (0:359)' * (2*pi/360);
    
    % Propagation vector at the angle of maximum energy
    ul = [sin(angles_rad( max_idx)); cos(angles_rad(max_idx))];
    
    % All possible propagation vectors
    ui_all = [sin(angles_rad), cos(angles_rad)]';
    
    % Dot products
    dot_uil = ul' * ui_all;
    
    % Weighted average of dot products
    planarity = sum(flux .* dot_uil') / sum(flux) .* 100;
end