function [ac_db, alpha, details] = full_racc_worst_case_contrast( ...
        w, HB, HD, etaB, etaD, deltaB, deltaD)
%FULL_RACC_WORST_CASE_CONTRAST Exact ATF-ball worst-case AC for a fixed w.
%
% For ||Delta H_Z||_F <= eta_Z, Delta H_Z*w can be any vector with
% Euclidean norm no larger than eta_Z*||w||_2. Consequently,
%
%   min ||(HB + Delta HB)w||_2^2 + deltaB*||w||_2^2
%       = max(||HB*w||_2 - etaB*||w||_2, 0)^2
%         + deltaB*||w||_2^2,
%
%   max ||(HD + Delta HD)w||_2^2 + deltaD*||w||_2^2
%       = (||HD*w||_2 + etaD*||w||_2)^2
%         + deltaD*||w||_2^2.
%
% This closed form is used only to score a recovered rank-one filter. It
% does not replace the lifted SDP used to design Full-RACC.

    if nargin < 6 || isempty(deltaB)
        deltaB = 0;
    end
    if nargin < 7 || isempty(deltaD)
        deltaD = 0;
    end

    validateattributes(w, {'numeric'}, {'column', 'nonempty'});
    validateattributes(etaB, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(etaD, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(deltaB, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(deltaD, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});

    MB = size(HB, 1);
    MD = size(HD, 1);
    if size(HB, 2) ~= numel(w) || size(HD, 2) ~= numel(w)
        error('FullRACC:DimensionMismatch', ...
            'HB and HD must have the same number of columns as w.');
    end

    w_norm = norm(w, 2);
    bright_nominal_norm = norm(HB * w, 2);
    dark_nominal_norm = norm(HD * w, 2);

    bright_worst_norm = max(bright_nominal_norm - etaB * w_norm, 0);
    dark_worst_norm = dark_nominal_norm + etaD * w_norm;
    bright_energy = bright_worst_norm^2 + deltaB * w_norm^2;
    dark_energy = dark_worst_norm^2 + deltaD * w_norm^2;

    numerator = MD * bright_energy;
    denominator = MB * dark_energy;
    if denominator == 0
        if numerator > 0
            alpha = Inf;
        else
            alpha = 0;
        end
    else
        alpha = numerator / denominator;
    end

    if alpha > 0
        ac_db = 10 * log10(alpha);
    elseif alpha == 0
        ac_db = -Inf;
    else
        ac_db = NaN;
    end

    details = struct( ...
        'bright_nominal_norm', bright_nominal_norm, ...
        'dark_nominal_norm', dark_nominal_norm, ...
        'bright_worst_energy', bright_energy, ...
        'dark_worst_energy', dark_energy, ...
        'deltaB', deltaB, ...
        'deltaD', deltaD, ...
        'w_norm', w_norm);
end
