function [filters, result] = vast_dft_endpoints( ...
        HB_ctrl, HD_ctrl, desired_BZ, para)
%VAST_DFT_ENDPOINTS Narrowband DFT-domain VAST endpoint filters.
%   FILTERS = VAST_DFT_ENDPOINTS(HB, HD, D, PARA) computes the VAST-NF
%   solutions independently at every frequency using the existing project
%   convention. The requested ranks are supplied in PARA.VAST.ranks.
%
%   For each frequency, the generalized eigenvectors U satisfy
%
%       U' * R_B * U = diag(lambda),  U' * R_D * U = I,
%
%   and the rank-V solution is
%
%       q_V = U_V * (U_V' * r_B ./ (mu + lambda_V)),
%
%   with r_B = H_B' * d. Rank V=1 is the ACC endpoint and V=L is the PM
%   endpoint of this narrowband VAST formulation.
%
%   Required fields:
%       para.vast.ranks  Integer vector in [1,L].
%
%   Optional fields:
%       para.vast.mu                     Default 1.
%       para.vast.dark_loading_ratio     Default 0 (paper SICER setting).
%       para.vast.fail_rcond_threshold   Default 0 (diagnostic only).
%
%   A positive dark_loading_ratio applies
%       R_D <- R_D + ratio*lambda_max(R_D)*I.
%   No automatic loading or fallback is applied.

    if nargin < 4 || ~isstruct(para) || ~isfield(para, 'vast') || ...
            ~isstruct(para.vast) || ~isfield(para.vast, 'ranks')
        error('VAST_DFT:MissingParameters', ...
            'Set para.vast.ranks explicitly.');
    end
    if ~isnumeric(HB_ctrl) || ~isnumeric(HD_ctrl) || ...
            ~isnumeric(desired_BZ) || any(~isfinite(HB_ctrl(:))) || ...
            any(~isfinite(HD_ctrl(:))) || any(~isfinite(desired_BZ(:)))
        error('VAST_DFT:InvalidInput', ...
            'ATFs and desired pressures must be finite numeric arrays.');
    end

    [bright_points, loudspeakers, num_frequencies] = size(HB_ctrl);
    [dark_points, dark_loudspeakers, dark_frequencies] = size(HD_ctrl);
    if bright_points < 1 || dark_points < 1 || ...
            dark_loudspeakers ~= loudspeakers || ...
            dark_frequencies ~= num_frequencies || ...
            ~isequal(size(desired_BZ), [bright_points, num_frequencies])
        error('VAST_DFT:DimensionMismatch', ...
            'ATF and desired-pressure dimensions are inconsistent.');
    end

    settings = local_settings(para.vast, loudspeakers);
    ranks = settings.ranks;
    num_ranks = numel(ranks);
    weights = complex(zeros(loudspeakers, num_frequencies, num_ranks));
    generalized_eigenvalues = nan(loudspeakers, num_frequencies);
    dark_rcond = nan(1, num_frequencies);
    applied_loading = nan(1, num_frequencies);
    endpoint_norms = nan(num_ranks, num_frequencies);
    full_rank_relative_residual = nan(1, num_frequencies);

    identity = eye(loudspeakers);
    for frequency_index = 1:num_frequencies
        HB = HB_ctrl(:, :, frequency_index);
        HD = HD_ctrl(:, :, frequency_index);
        desired = desired_BZ(:, frequency_index);
        R_B = local_hermitian(HB' * HB);
        R_D = local_hermitian(HD' * HD);
        r_B = HB' * desired;

        if settings.dark_loading_ratio > 0
            dark_scale = max(real(eig(R_D, 'vector')));
            loading = settings.dark_loading_ratio * dark_scale;
        else
            loading = 0;
        end
        R_D_loaded = local_hermitian(R_D + loading * identity);
        dark_rcond(frequency_index) = rcond(R_D_loaded);
        applied_loading(frequency_index) = loading;
        if dark_rcond(frequency_index) < settings.fail_rcond_threshold
            error('VAST_DFT:DarkMatrixCondition', ...
                ['R_D has rcond %.3e at frequency index %d, below the ', ...
                 'declared threshold %.3e.'], dark_rcond(frequency_index), ...
                frequency_index, settings.fail_rcond_threshold);
        end

        [U, lambda] = jdiag(R_B, R_D_loaded, 'vector');
        lambda = real(lambda(:));
        generalized_eigenvalues(:, frequency_index) = lambda;
        projected_target = U' * r_B;

        for rank_index = 1:num_ranks
            rank_value = ranks(rank_index);
            U_part = U(:, 1:rank_value);
            coefficients = projected_target(1:rank_value) ./ ...
                (settings.mu + lambda(1:rank_value));
            q = U_part * coefficients;
            if any(~isfinite(q)) || norm(q) <= 0
                error('VAST_DFT:InvalidFilter', ...
                    'Invalid rank-%d filter at frequency index %d.', ...
                    rank_value, frequency_index);
            end
            weights(:, frequency_index, rank_index) = q;
            endpoint_norms(rank_index, frequency_index) = norm(q);
        end

        full_rank_index = find(ranks == loudspeakers, 1);
        if ~isempty(full_rank_index)
            q_full = weights(:, frequency_index, full_rank_index);
            normal_matrix = R_B + settings.mu * R_D_loaded;
            residual_scale = max(norm(r_B), eps);
            full_rank_relative_residual(frequency_index) = ...
                norm(normal_matrix * q_full - r_B) / residual_scale;
        end
    end

    filters = struct();
    for rank_index = 1:num_ranks
        field_name = sprintf('V%d', ranks(rank_index));
        filters.(field_name) = weights(:, :, rank_index);
    end

    result = struct();
    result.algorithm = 'DFT-domain VAST-NF';
    result.ranks = ranks;
    result.mu = settings.mu;
    result.dark_loading_ratio = settings.dark_loading_ratio;
    result.applied_dark_loading = applied_loading;
    result.dark_rcond = dark_rcond;
    result.generalized_eigenvalues = generalized_eigenvalues;
    result.filter_norms = endpoint_norms;
    result.full_rank_relative_residual = full_rank_relative_residual;
    result.endpoint_labels = local_labels(ranks, loudspeakers);
end

function settings = local_settings(source, loudspeakers)
    settings.ranks = source.ranks(:).';
    validateattributes(settings.ranks, {'numeric'}, ...
        {'vector', 'integer', 'positive', '<=', loudspeakers});
    if numel(unique(settings.ranks)) ~= numel(settings.ranks)
        error('VAST_DFT:DuplicateRanks', 'VAST ranks must be unique.');
    end
    settings.mu = local_field(source, 'mu', 1);
    settings.dark_loading_ratio = local_field( ...
        source, 'dark_loading_ratio', 0);
    settings.fail_rcond_threshold = local_field( ...
        source, 'fail_rcond_threshold', 0);
    validateattributes(settings.mu, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(settings.dark_loading_ratio, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(settings.fail_rcond_threshold, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
end

function value = local_field(source, name, default_value)
    if isfield(source, name)
        value = source.(name);
    else
        value = default_value;
    end
end

function matrix = local_hermitian(matrix)
    matrix = (matrix + matrix') / 2;
end

function labels = local_labels(ranks, loudspeakers)
    labels = strings(size(ranks));
    for k = 1:numel(ranks)
        if ranks(k) == 1
            labels(k) = "ACC endpoint";
        elseif ranks(k) == loudspeakers
            labels(k) = "PM endpoint";
        else
            labels(k) = sprintf('V=%d', ranks(k));
        end
    end
end
