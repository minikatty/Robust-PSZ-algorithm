function res = NoCrossTermWCRACC(HB_ctrl, HD_ctrl, para, fre_indices)
%NOCROSSTERMWCRACC WCRACC after removing the nominal-error cross terms.
%
% The ablated correlation model is
%
%   Rtilde_Z = R_Z + DeltaH_Z^H*DeltaH_Z,
%   ||DeltaH_Z||_F <= eta_Z.
%
% Since DeltaH_Z^H*DeltaH_Z is positive semidefinite, the worst bright-zone
% energy is attained at DeltaH_B=0, while the worst dark-zone energy adds
% eta_D^2*||w||_2^2. The resulting generalized eigenproblem is therefore
%
%   max_w (w^H (R_B+delta_B I) w) /
%         (w^H (R_D+(delta_D+eta_D^2) I) w).
%
% This strict NoCT model is the direct ablation counterpart of FullWCRACC.
% It is intentionally separate from the historical wcACC.m implementation,
% which uses an unstructured correlation-matrix uncertainty ball.

    if ndims(HB_ctrl) == 2
        HB_ctrl = reshape(HB_ctrl, size(HB_ctrl, 1), size(HB_ctrl, 2), 1);
    end
    if ndims(HD_ctrl) == 2
        HD_ctrl = reshape(HD_ctrl, size(HD_ctrl, 1), size(HD_ctrl, 2), 1);
    end

    [MB, L, num_freqs] = size(HB_ctrl);
    [MD, LD, num_freqs_D] = size(HD_ctrl);
    if LD ~= L || num_freqs_D ~= num_freqs
        error('NoCTWCRACC:DimensionMismatch', ...
            ['HB_ctrl and HD_ctrl must have compatible loudspeaker ', ...
             'and frequency dimensions.']);
    end
    if nargin < 4 || isempty(fre_indices)
        fre_indices = 1:num_freqs;
    end
    if any(fre_indices < 1) || any(fre_indices > num_freqs) || ...
            any(fre_indices ~= round(fre_indices))
        error('NoCTWCRACC:InvalidFrequencyIndices', ...
            'fre_indices must contain valid integer frequency indices.');
    end

    loading_ratio = local_loading_ratio(para);
    w_all = zeros(L, num_freqs);
    etaB_all = nan(1, num_freqs);
    etaD_all = nan(1, num_freqs);
    deltaB_all = nan(1, num_freqs);
    deltaD_all = nan(1, num_freqs);
    model_alpha = nan(1, num_freqs);
    model_ac_db = nan(1, num_freqs);
    full_ball_alpha = nan(1, num_freqs);
    full_ball_ac_db = nan(1, num_freqs);

    for i = fre_indices(:).'
        HB = HB_ctrl(:, :, i);
        HD = HD_ctrl(:, :, i);
        etaB = local_resolve_eta(para, 'B', i, HB);
        etaD = local_resolve_eta(para, 'D', i, HD);
        deltaB = loading_ratio * norm(HB, 2)^2;
        deltaD = loading_ratio * norm(HD, 2)^2;

        data_scale = max(norm([HB; HD], 'fro'), realmin);
        HBs = HB / data_scale;
        HDs = HD / data_scale;
        etaDs = etaD / data_scale;
        deltaBs = deltaB / data_scale^2;
        deltaDs = deltaD / data_scale^2;

        RB = HBs' * HBs;
        RD = HDs' * HDs;
        RB = (RB + RB') / 2;
        RD = (RD + RD') / 2;
        A = RB + deltaBs * eye(L);
        B = RD + (deltaDs + etaDs^2) * eye(L);

        numerical_floor = max(1e-12 * max(real(trace(B)) / L, 1), eps);
        [V, D] = eig(A, B + numerical_floor * eye(L));
        [~, max_index] = max(real(diag(D)));
        w = V(:, max_index);
        w = w / max(norm(w), realmin);

        numerator = real(w' * A * w);
        denominator = real(w' * B * w);
        alpha = (MD / MB) * numerator / max(denominator, realmin);
        [full_db, full_alpha] = full_racc_worst_case_contrast( ...
            w, HB, HD, etaB, etaD, deltaB, deltaD);

        w_all(:, i) = w;
        etaB_all(i) = etaB;
        etaD_all(i) = etaD;
        deltaB_all(i) = deltaB;
        deltaD_all(i) = deltaD;
        model_alpha(i) = alpha;
        model_ac_db(i) = 10 * log10(max(alpha, realmin));
        full_ball_alpha(i) = full_alpha;
        full_ball_ac_db(i) = full_db;
    end

    res = struct();
    res.w = w_all;
    res.scale = 1;
    res.algorithm = 'NoCT-WCRACC';
    res.eta = struct('B', etaB_all, 'D', etaD_all);
    res.gamma = struct('B', etaB_all.^2, 'D', etaD_all.^2);
    res.diagonal_loading = struct('B', deltaB_all, 'D', deltaD_all, ...
        'ratio', loading_ratio);
    res.model_alpha = model_alpha;
    res.model_ac_db = model_ac_db;
    res.full_atf_worst_case_alpha = full_ball_alpha;
    res.full_atf_worst_case_ac_db = full_ball_ac_db;
end

function ratio = local_loading_ratio(para)
    ratio = 0;
    if isfield(para, 'full_racc') && ...
            isfield(para.full_racc, 'diagonal_loading_ratio')
        ratio = para.full_racc.diagonal_loading_ratio;
    end
    validateattributes(ratio, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
end

function eta = local_resolve_eta(para, zone, freq_index, H)
    if isfield(para, 'full_racc') && isfield(para.full_racc, 'eta') && ...
            isfield(para.full_racc.eta, zone)
        source = para.full_racc.eta.(zone);
        eta = local_pick_value(source, freq_index);
    elseif isfield(para, 'epsilon') && isfield(para.epsilon, zone)
        eta = local_pick_value(para.epsilon.(zone), freq_index);
    elseif isfield(para, 'full_racc') && ...
            isfield(para.full_racc, 'relative_eta')
        source = para.full_racc.relative_eta;
        if isstruct(source)
            source = source.(zone);
        end
        eta = local_pick_value(source, freq_index) * norm(H, 'fro');
    else
        error('NoCTWCRACC:MissingUncertaintyRadius', ...
            ['Provide para.full_racc.eta.%s, para.epsilon.%s, or ', ...
             'para.full_racc.relative_eta.'], zone, zone);
    end
    validateattributes(eta, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
end

function value = local_pick_value(source, frequency_index)
    if isscalar(source)
        value = source;
    elseif numel(source) >= frequency_index
        value = source(frequency_index);
    else
        error('NoCTWCRACC:UncertaintyLengthMismatch', ...
            'No uncertainty radius is available for frequency index %d.', ...
            frequency_index);
    end
end
