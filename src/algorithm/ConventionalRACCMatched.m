function res = ConventionalRACCMatched(HB_ctrl, HD_ctrl, para, fre_indices)
%CONVENTIONALRACCMATCHED Correlation-level RACC with gamma_Z = eta_Z^2.
%
% This is a parameter-matched ablation baseline for FullCrossTermRACC. The
% ATF uncertainty radii eta_B and eta_D are resolved using the same rules as
% FullCrossTermRACC, then mapped to correlation-level radii by
% gamma_Z = eta_Z^2. This function deliberately does not modify wcACC.m,
% whose historical 0.0075*||R||_F setting is used by existing experiments.
% An optional common regularization is set by
% para.full_racc.diagonal_loading_ratio, with
% delta_Z=tau*lambda_max(R_Z).

    if ndims(HB_ctrl) == 2
        HB_ctrl = reshape(HB_ctrl, size(HB_ctrl, 1), size(HB_ctrl, 2), 1);
    end
    if ndims(HD_ctrl) == 2
        HD_ctrl = reshape(HD_ctrl, size(HD_ctrl, 1), size(HD_ctrl, 2), 1);
    end

    [~, L, num_freqs] = size(HB_ctrl);
    if size(HD_ctrl, 2) ~= L || size(HD_ctrl, 3) ~= num_freqs
        error('ConventionalRACCMatched:DimensionMismatch', ...
            'HB_ctrl and HD_ctrl must have compatible loudspeaker/frequency dimensions.');
    end
    if nargin < 4 || isempty(fre_indices)
        fre_indices = 1:num_freqs;
    end

    w_all = zeros(L, num_freqs);
    etaB_all = nan(1, num_freqs);
    etaD_all = nan(1, num_freqs);
    gammaB_all = nan(1, num_freqs);
    gammaD_all = nan(1, num_freqs);
    model_alpha = nan(1, num_freqs);
    model_ac_db = nan(1, num_freqs);
    full_alpha = nan(1, num_freqs);
    full_ac_db = nan(1, num_freqs);
    deltaB_all = nan(1, num_freqs);
    deltaD_all = nan(1, num_freqs);
    loading_ratio = local_loading_ratio(para);

    for i = fre_indices(:).'
        HB = HB_ctrl(:, :, i);
        HD = HD_ctrl(:, :, i);
        etaB = local_resolve_eta(para, 'B', i, HB);
        etaD = local_resolve_eta(para, 'D', i, HD);
        deltaB = loading_ratio * norm(HB, 2)^2;
        deltaD = loading_ratio * norm(HD, 2)^2;

        % Scaling H and eta together leaves the generalized eigenvector
        % unchanged and avoids poorly scaled correlation matrices.
        data_scale = max(norm([HB; HD], 'fro'), 1e-12);
        HBs = HB / data_scale;
        HDs = HD / data_scale;
        etaBs = etaB / data_scale;
        etaDs = etaD / data_scale;
        deltaBs = deltaB / data_scale^2;
        deltaDs = deltaD / data_scale^2;

        RB = (HBs' * HBs + (HBs' * HBs)') / 2;
        RD = (HDs' * HDs + (HDs' * HDs)') / 2;
        gammaB = etaBs^2;
        gammaD = etaDs^2;
        A = RB + (deltaBs - gammaB) * eye(L);
        B = RD + (deltaDs + gammaD) * eye(L);
        A = (A + A') / 2;
        B = (B + B') / 2;

        % B is positive definite when etaD > 0. Add only a numerical floor
        % for the zero-radius/singular-RD edge case.
        b_floor = max(1e-12 * max(real(trace(B)) / L, 1), eps);
        B = B + b_floor * eye(L);
        [V, D] = eig(A, B);
        eig_values = real(diag(D));
        [~, max_idx] = max(eig_values);
        w = V(:, max_idx);
        if norm(w) > 0
            w = w / norm(w);
        end

        MB = size(HB, 1);
        MD = size(HD, 1);
        corr_num = real(w' * A * w);
        corr_den = real(w' * B * w);
        alpha_corr = (MD / MB) * corr_num / max(corr_den, realmin);
        [ac_full_db_i, alpha_full_i] = ...
            full_racc_worst_case_contrast( ...
            w, HB, HD, etaB, etaD, deltaB, deltaD);

        w_all(:, i) = w;
        etaB_all(i) = etaB;
        etaD_all(i) = etaD;
        gammaB_all(i) = etaB^2;
        gammaD_all(i) = etaD^2;
        model_alpha(i) = alpha_corr;
        model_ac_db(i) = 10 * log10(max(alpha_corr, realmin));
        full_alpha(i) = alpha_full_i;
        full_ac_db(i) = ac_full_db_i;
        deltaB_all(i) = deltaB;
        deltaD_all(i) = deltaD;
    end

    res = struct();
    res.w = w_all;
    res.scale = 1;
    res.eta = struct('B', etaB_all, 'D', etaD_all);
    res.gamma = struct('B', gammaB_all, 'D', gammaD_all);
    res.diagonal_loading = struct('B', deltaB_all, 'D', deltaD_all, ...
        'ratio', loading_ratio);
    res.correlation_model_alpha = model_alpha;
    res.correlation_model_ac_db = model_ac_db;
    res.full_atf_worst_case_alpha = full_alpha;
    res.full_atf_worst_case_ac_db = full_ac_db;
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

function eta = local_resolve_eta(para, zone, freq_idx, H)
    if isfield(para, 'full_racc') && isfield(para.full_racc, 'eta') && ...
            isfield(para.full_racc.eta, zone)
        eta_source = para.full_racc.eta.(zone);
    elseif isfield(para, 'epsilon') && isfield(para.epsilon, zone)
        eta_source = para.epsilon.(zone);
    elseif isfield(para, 'full_racc') && isfield(para.full_racc, 'relative_eta')
        relative_source = para.full_racc.relative_eta;
        if isstruct(relative_source)
            relative_source = relative_source.(zone);
        end
        eta = local_pick_value(relative_source, freq_idx) * norm(H, 'fro');
        local_validate_eta(eta, zone, freq_idx);
        return;
    else
        error('FullRACC:MissingUncertaintyRadius', ...
            ['Provide para.full_racc.eta.%s, para.epsilon.%s, or ', ...
             'para.full_racc.relative_eta.'], zone, zone);
    end
    eta = local_pick_value(eta_source, freq_idx);
    local_validate_eta(eta, zone, freq_idx);
end

function value = local_pick_value(source, freq_idx)
    if isscalar(source)
        value = source;
    elseif numel(source) >= freq_idx
        value = source(freq_idx);
    else
        error('FullRACC:UncertaintyLengthMismatch', ...
            'The uncertainty-radius vector does not contain frequency index %d.', freq_idx);
    end
end

function local_validate_eta(eta, zone, freq_idx)
    if ~isnumeric(eta) || ~isscalar(eta) || ~isreal(eta) || ...
            ~isfinite(eta) || eta < 0
        error('FullRACC:InvalidUncertaintyRadius', ...
            'eta_%s at frequency index %d must be a finite nonnegative scalar.', ...
            zone, freq_idx);
    end
end
