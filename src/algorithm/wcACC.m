function res = wcACC(HB_ctrl, HD_ctrl, wcACC_para)
%WCACC Conventional correlation-level worst-case robust ACC.
%
% When an ATF uncertainty radius is supplied, this implementation uses the
% common mapping gamma_Z=eta_Z^2. Historical scripts that do not provide an
% explicit/relative ATF radius retain the former 0.0075*||R_B||_F setting.

    frePoint = wcACC_para.target_freqs;
    num_frePoint = length(frePoint);
    L = size(HB_ctrl, 2);
    w = zeros(L, num_frePoint);
    gammaB_all = nan(1, num_frePoint);
    gammaD_all = nan(1, num_frePoint);
    etaB_all = nan(1, num_frePoint);
    etaD_all = nan(1, num_frePoint);
    deltaB_all = nan(1, num_frePoint);
    deltaD_all = nan(1, num_frePoint);
    loading_ratio = local_loading_ratio(wcACC_para);
    use_common_radius = local_has_common_radius(wcACC_para);

    for i = 1:num_frePoint
        HB = squeeze(HB_ctrl(:, :, i));
        HD = squeeze(HD_ctrl(:, :, i));
        RB = HB' * HB;
        RD = HD' * HD;
        RB = (RB + RB') / 2;
        RD = (RD + RD') / 2;

        if use_common_radius
            etaB = local_resolve_eta(wcACC_para, 'B', i, HB);
            etaD = local_resolve_eta(wcACC_para, 'D', i, HD);
            GammaB = etaB^2;
            GammaD = etaD^2;
        else
            GammaB = 0.75 * norm(RB, 'fro') / 100;
            GammaD = GammaB;
            etaB = sqrt(GammaB);
            etaD = sqrt(GammaD);
        end

        deltaB = loading_ratio * norm(HB, 2)^2;
        deltaD = loading_ratio * norm(HD, 2)^2;
        numerator_matrix = RB + (deltaB - GammaB) * eye(L);
        denominator_matrix = RD + (deltaD + GammaD) * eye(L);
        denominator_matrix = (denominator_matrix + denominator_matrix') / 2;
        numerical_floor = max(1e-12 * ...
            max(real(trace(denominator_matrix)) / L, 1), eps);

        A = (denominator_matrix + numerical_floor * eye(L)) \ ...
            ((numerator_matrix + numerator_matrix') / 2);
        [~, w(:, i)] = MaxEigenvector(A);

        gammaB_all(i) = GammaB;
        gammaD_all(i) = GammaD;
        etaB_all(i) = etaB;
        etaD_all(i) = etaD;
        deltaB_all(i) = deltaB;
        deltaD_all(i) = deltaD;
    end

    res.w = w;
    res.scale = wcACC_para.scale;
    res.eta = struct('B', etaB_all, 'D', etaD_all);
    res.gamma = struct('B', gammaB_all, 'D', gammaD_all);
    res.diagonal_loading = struct('B', deltaB_all, 'D', deltaD_all, ...
        'ratio', loading_ratio);
    res.uses_common_atf_radius = use_common_radius;
end

function tf = local_has_common_radius(para)
    tf = (isfield(para, 'full_racc') && ...
        (isfield(para.full_racc, 'eta') || ...
         isfield(para.full_racc, 'relative_eta'))) || ...
        (isfield(para, 'epsilon') && isfield(para.epsilon, 'B') && ...
         isfield(para.epsilon, 'D'));
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

function eta = local_resolve_eta(para, zone, frequency_index, H)
    if isfield(para, 'full_racc') && isfield(para.full_racc, 'eta') && ...
            isfield(para.full_racc.eta, zone)
        eta = local_pick_value(para.full_racc.eta.(zone), frequency_index);
    elseif isfield(para, 'epsilon') && isfield(para.epsilon, zone)
        eta = local_pick_value(para.epsilon.(zone), frequency_index);
    elseif isfield(para, 'full_racc') && ...
            isfield(para.full_racc, 'relative_eta')
        source = para.full_racc.relative_eta;
        if isstruct(source)
            source = source.(zone);
        end
        eta = local_pick_value(source, frequency_index) * norm(H, 'fro');
    else
        error('wcACC:MissingUncertaintyRadius', ...
            'The common ATF uncertainty radius is incomplete.');
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
        error('wcACC:UncertaintyLengthMismatch', ...
            'No uncertainty radius is available for frequency index %d.', ...
            frequency_index);
    end
end
