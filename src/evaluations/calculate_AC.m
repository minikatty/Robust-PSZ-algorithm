function AC = calculate_AC(w, H_B, H_D)
% CALCULATE_AC - Computes the Acoustic Contrast (AC) in dB.
% The contrast is the ratio of mean bright-zone to mean dark-zone energy.
% Do not add a fixed epsilon to either physical energy: doing so breaks the
% required invariance AC(c*w)=AC(w) and biases algorithms whose filters have
% different overall scales.
    M_B = size(H_B, 1);
    M_D = size(H_D, 1);

    bright_energy = norm(H_B * w, 2)^2;
    dark_energy = norm(H_D * w, 2)^2;
    mu = (M_D * bright_energy) / ...
        max(M_B * dark_energy, realmin);
    AC = 10 * log10(max(real(mu), realmin));
end
