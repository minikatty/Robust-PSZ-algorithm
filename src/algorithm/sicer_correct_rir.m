function [corrected_rir, info] = sicer_correct_rir( ...
        old_rir, old_sound_speed, new_sound_speed, fs, options)
%SICER_CORRECT_RIR Correct RIRs for a uniform sound-speed change.
%   CORRECTED_RIR = SICER_CORRECT_RIR(OLD_RIR, C_OLD, C_NEW, FS)
%   applies the SICER model
%
%       h_new(m) = (1/beta) * sum_n h_old(n) sinc(m/beta - n),
%       beta = C_OLD / C_NEW,
%
%   along the last dimension of OLD_RIR. The output has the same size as
%   the input. For C_NEW > C_OLD (beta < 1), the old RIR is first low-pass
%   filtered at beta*FS/2 using the paper's order-100 FIR convention and
%   the linear-phase delay is compensated.
%
%   OPTIONS fields:
%       lowpass_order       Even FIR order (default 100).
%       apply_antialias     Logical scalar (default true).
%       chunk_size          Number of paths per matrix product (default 512).
%       beta_identity_tol   Identity tolerance (default 1e-12).
%
%   The implementation forms the interpolation matrix once per beta and
%   applies it to path chunks. It is intended as the reference time-domain
%   implementation against which faster frequency-domain realizations can
%   be checked.

    if nargin < 5 || isempty(options)
        options = struct();
    end
    validateattributes(old_rir, {'numeric'}, {'nonempty', 'finite'});
    validateattributes(old_sound_speed, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(new_sound_speed, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(fs, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});

    settings = local_settings(options);
    beta = double(old_sound_speed) / double(new_sound_speed);
    original_size = size(old_rir);
    num_samples = original_size(end);
    num_paths = numel(old_rir) / num_samples;

    path_matrix = reshape(permute(old_rir, ...
        [numel(original_size), 1:numel(original_size)-1]), ...
        num_samples, num_paths);

    total_timer = tic;
    lowpass_applied = settings.apply_antialias && beta < 1 - ...
        settings.beta_identity_tol;
    if lowpass_applied
        normalized_cutoff = beta;
        lowpass_coefficients = fir1(settings.lowpass_order, ...
            normalized_cutoff, 'low');
        filter_delay = settings.lowpass_order / 2;
    else
        normalized_cutoff = NaN;
        lowpass_coefficients = 1;
        filter_delay = 0;
    end

    if abs(beta - 1) <= settings.beta_identity_tol
        corrected_matrix = path_matrix;
        interpolation_seconds = 0;
        kernel_bytes = 0;
    else
        old_indices = (0:num_samples-1).';
        new_positions = (0:num_samples-1) / beta;
        interpolation_arguments = new_positions - old_indices;
        interpolation_kernel = local_normalized_sinc( ...
            interpolation_arguments);
        kernel_info = whos('interpolation_kernel');
        kernel_bytes = kernel_info.bytes;

        corrected_matrix = zeros(size(path_matrix), 'like', path_matrix);
        interpolation_timer = tic;
        for first_path = 1:settings.chunk_size:num_paths
            last_path = min(first_path + settings.chunk_size - 1, ...
                num_paths);
            current = path_matrix(:, first_path:last_path);
            if lowpass_applied
                current = filter(lowpass_coefficients, 1, current, [], 1);
                if filter_delay > 0
                    current = [current(filter_delay+1:end, :); ...
                        zeros(filter_delay, size(current, 2), ...
                        'like', current)];
                end
            end
            corrected_matrix(:, first_path:last_path) = ...
                (interpolation_kernel.' * current) / beta;
        end
        interpolation_seconds = toc(interpolation_timer);
    end

    permuted_size = [num_samples, original_size(1:end-1)];
    corrected_permuted = reshape(corrected_matrix, permuted_size);
    corrected_rir = ipermute(corrected_permuted, ...
        [numel(original_size), 1:numel(original_size)-1]);

    info = struct();
    info.method = 'SICER time-domain sinc interpolation';
    info.old_sound_speed_mps = double(old_sound_speed);
    info.new_sound_speed_mps = double(new_sound_speed);
    info.beta = beta;
    info.num_samples = num_samples;
    info.num_paths = num_paths;
    info.output_size = size(corrected_rir);
    info.antialias_lowpass_applied = lowpass_applied;
    info.lowpass_order = settings.lowpass_order;
    info.normalized_cutoff = normalized_cutoff;
    info.filter_delay_samples = filter_delay;
    info.chunk_size = settings.chunk_size;
    info.interpolation_kernel_bytes = kernel_bytes;
    info.interpolation_seconds = interpolation_seconds;
    info.total_seconds = toc(total_timer);
end

function settings = local_settings(options)
    settings.lowpass_order = local_field(options, 'lowpass_order', 100);
    settings.apply_antialias = local_field(options, 'apply_antialias', true);
    settings.chunk_size = local_field(options, 'chunk_size', 512);
    settings.beta_identity_tol = local_field( ...
        options, 'beta_identity_tol', 1e-12);

    validateattributes(settings.lowpass_order, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative', 'even'});
    validateattributes(settings.apply_antialias, {'logical', 'numeric'}, ...
        {'scalar'});
    settings.apply_antialias = logical(settings.apply_antialias);
    validateattributes(settings.chunk_size, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(settings.beta_identity_tol, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
end

function value = local_field(source, name, default_value)
    if isfield(source, name)
        value = source.(name);
    else
        value = default_value;
    end
end

function values = local_normalized_sinc(arguments)
    values = ones(size(arguments));
    nonzero = arguments ~= 0;
    scaled = pi * arguments(nonzero);
    values(nonzero) = sin(scaled) ./ scaled;
end
