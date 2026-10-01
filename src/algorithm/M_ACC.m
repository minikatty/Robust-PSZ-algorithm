function [w, result] = M_ACC(HB_ctrl, HD_ctrl, para)
%M_ACC Modified acoustic contrast control with bright-zone uniformity.
%   [W, RESULT] = M_ACC(HB_CTRL, HD_CTRL, PARA) solves, independently at
%   every frequency,
%
%       maximize_q  (q' * R_B * q) / ...
%           (q' * (R_D + alpha * R_SUC + beta * I) * q),
%
%   where R_B = H_B' * H_B, R_D = H_D' * H_D, and
%
%       H_SUC = H_B - ones(M_B, 1) * mean(H_B, 1),
%       R_SUC = H_SUC' * H_SUC.
%
%   This is the frequency-domain modified ACC formulation in Li et al.,
%   Applied Acoustics 253 (2026) 111374. The returned column at each
%   frequency has unit Euclidean norm, consistent with the project's ACC
%   implementations. Its arbitrary global phase is fixed deterministically.
%
%   Required parameter:
%       para.macc.alpha       Scalar or one value per frequency.
%
%   Optional parameters:
%       para.macc.beta        Scalar or one value per frequency (default
%                             1e-14, as reported in the paper).
%       para.macc.beta_mode   'absolute' (paper setting, default) or
%                             'relative_dark_maxeig'. In the relative mode,
%                             the applied loading is beta*lambda_max(R_D).
%
%   M-ACC is a nominal spatial-uniformity method. It is not an uncertainty-
%   robust design, and alpha must be selected without held-out evaluation
%   ATFs.

    if nargin < 3 || ~isstruct(para)
        error('M_ACC:MissingParameters', ...
            'A parameter structure containing para.macc.alpha is required.');
    end
    if ~isnumeric(HB_ctrl) || ~isnumeric(HD_ctrl) || ...
            any(~isfinite(HB_ctrl(:))) || any(~isfinite(HD_ctrl(:)))
        error('M_ACC:InvalidATF', ...
            'Bright- and dark-zone ATFs must be finite numeric arrays.');
    end

    [bright_points, loudspeakers, num_frequencies] = size(HB_ctrl);
    [dark_points, dark_loudspeakers, dark_frequencies] = size(HD_ctrl);
    if bright_points < 1 || dark_points < 1 || loudspeakers < 1 || ...
            dark_loudspeakers ~= loudspeakers || ...
            dark_frequencies ~= num_frequencies
        error('M_ACC:DimensionMismatch', ...
            ['HB_ctrl and HD_ctrl must contain at least one point and ', ...
             'share their loudspeaker and frequency dimensions.']);
    end

    settings = local_read_settings(para, num_frequencies);
    alpha = settings.alpha;
    beta = settings.beta;

    w = complex(zeros(loudspeakers, num_frequencies));
    applied_beta = nan(1, num_frequencies);
    objective_value = nan(1, num_frequencies);
    bright_energy = nan(1, num_frequencies);
    dark_energy = nan(1, num_frequencies);
    uniformity_energy = nan(1, num_frequencies);
    denominator_energy = nan(1, num_frequencies);
    denominator_rcond = nan(1, num_frequencies);
    stationarity_residual = nan(1, num_frequencies);

    identity = eye(loudspeakers);
    for frequency_index = 1:num_frequencies
        HB = HB_ctrl(:, :, frequency_index);
        HD = HD_ctrl(:, :, frequency_index);

        mean_HB = ones(bright_points, 1) * mean(HB, 1);
        H_SUC = HB - mean_HB;

        R_B = local_hermitian(HB' * HB);
        R_D = local_hermitian(HD' * HD);
        R_SUC = local_hermitian(H_SUC' * H_SUC);

        switch settings.beta_mode
            case 'absolute'
                beta_at_frequency = beta(frequency_index);
            case 'relative_dark_maxeig'
                dark_eigenvalues = eig(R_D, 'vector');
                dark_scale = max(real(dark_eigenvalues));
                beta_at_frequency = beta(frequency_index) * dark_scale;
            otherwise
                error('M_ACC:InternalBetaMode', ...
                    'Unsupported beta mode: %s.', settings.beta_mode);
        end

        Q = local_hermitian(R_D + ...
            alpha(frequency_index) * R_SUC + ...
            beta_at_frequency * identity);

        % Q must be positive definite for the stated regularized quotient.
        % Cholesky whitening preserves the Hermitian eigenproblem and avoids
        % explicitly forming inv(Q) or Q\R_B.
        [lower_factor, chol_flag] = chol(Q, 'lower');
        if chol_flag ~= 0
            error('M_ACC:NonPositiveDefiniteDenominator', ...
                ['The M-ACC denominator is not positive definite at ', ...
                 'frequency index %d. Increase beta or use ', ...
                 'relative_dark_maxeig loading.'], frequency_index);
        end

        whitened_bright = lower_factor \ R_B / lower_factor';
        whitened_bright = local_hermitian(whitened_bright);
        [eigenvectors, eigenvalues] = eig(whitened_bright, 'vector');
        [~, dominant_index] = max(real(eigenvalues));
        q = lower_factor' \ eigenvectors(:, dominant_index);

        q_norm = norm(q, 2);
        if ~isfinite(q_norm) || q_norm <= 0
            error('M_ACC:InvalidEigenvector', ...
                'Failed to recover a finite M-ACC vector at index %d.', ...
                frequency_index);
        end
        q = q / q_norm;
        q = local_fix_global_phase(q);
        w(:, frequency_index) = q;

        bright_value = max(real(q' * R_B * q), 0);
        dark_value = max(real(q' * R_D * q), 0);
        uniformity_value = max(real(q' * R_SUC * q), 0);
        denominator_value = real(q' * Q * q);
        quotient_value = bright_value / denominator_value;
        residual_denominator = max(norm(R_B * q, 2), eps);
        residual_value = norm(R_B * q - quotient_value * Q * q, 2) / ...
            residual_denominator;

        applied_beta(frequency_index) = beta_at_frequency;
        objective_value(frequency_index) = quotient_value;
        bright_energy(frequency_index) = bright_value;
        dark_energy(frequency_index) = dark_value;
        uniformity_energy(frequency_index) = uniformity_value;
        denominator_energy(frequency_index) = denominator_value;
        denominator_rcond(frequency_index) = rcond(Q);
        stationarity_residual(frequency_index) = residual_value;
    end

    result = struct();
    result.algorithm = 'M_ACC';
    result.alpha = alpha;
    result.beta_requested = beta;
    result.beta_mode = settings.beta_mode;
    result.beta_applied = applied_beta;
    result.objective_value = objective_value;
    result.bright_energy = bright_energy;
    result.dark_energy = dark_energy;
    result.uniformity_energy = uniformity_energy;
    result.denominator_energy = denominator_energy;
    result.denominator_rcond = denominator_rcond;
    result.stationarity_residual = stationarity_residual;
    result.normalization = 'unit_2_norm';
end

function settings = local_read_settings(para, num_frequencies)
    if ~isfield(para, 'macc') || ~isstruct(para.macc) || ...
            ~isfield(para.macc, 'alpha')
        error('M_ACC:MissingAlpha', ...
            ['Set para.macc.alpha explicitly. It must be chosen from ', ...
             'nominal control data, not held-out evaluation ATFs.']);
    end

    settings.alpha = local_expand_parameter( ...
        para.macc.alpha, num_frequencies, 'alpha');
    if isfield(para.macc, 'beta')
        beta_source = para.macc.beta;
    else
        beta_source = 1e-14;
    end
    settings.beta = local_expand_parameter( ...
        beta_source, num_frequencies, 'beta');

    if isfield(para.macc, 'beta_mode')
        beta_mode = para.macc.beta_mode;
        if ~(ischar(beta_mode) || (isstring(beta_mode) && isscalar(beta_mode)))
            error('M_ACC:InvalidBetaMode', ...
                'para.macc.beta_mode must be a character vector or string scalar.');
        end
        settings.beta_mode = lower(strtrim(char(beta_mode)));
    else
        settings.beta_mode = 'absolute';
    end
    valid_modes = {'absolute', 'relative_dark_maxeig'};
    if ~any(strcmp(settings.beta_mode, valid_modes))
        error('M_ACC:InvalidBetaMode', ...
            ['para.macc.beta_mode must be ''absolute'' or ', ...
             '''relative_dark_maxeig''.']);
    end
end

function values = local_expand_parameter(source, num_frequencies, name)
    validateattributes(source, {'numeric'}, ...
        {'vector', 'real', 'finite', 'nonnegative'}, mfilename, name);
    if isscalar(source)
        values = repmat(double(source), 1, num_frequencies);
    elseif numel(source) == num_frequencies
        values = reshape(double(source), 1, num_frequencies);
    else
        error('M_ACC:ParameterLengthMismatch', ...
            '%s must be scalar or contain one value per frequency.', name);
    end
end

function matrix = local_hermitian(matrix)
    matrix = (matrix + matrix') / 2;
end

function vector = local_fix_global_phase(vector)
    [~, reference_index] = max(abs(vector));
    reference_value = vector(reference_index);
    if abs(reference_value) > 0
        vector = vector * exp(-1i * angle(reference_value));
    end
end
