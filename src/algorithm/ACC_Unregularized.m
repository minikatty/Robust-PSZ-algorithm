function w = ACC_Unregularized(HB_ctrl, HD_ctrl)
%ACC_UNREGULARIZED Conventional ACC without diagonal loading.
%   This baseline deliberately solves R_D \ R_B directly. At frequencies
%   where R_D is ill-conditioned, MATLAB's singular/nearly-singular warning
%   is expected and is suppressed locally. No diagonal loading, pseudoinverse,
%   or other regularization is introduced.

    [~, loudspeakers, num_frequencies] = size(HB_ctrl);
    if size(HD_ctrl, 2) ~= loudspeakers || ...
            size(HD_ctrl, 3) ~= num_frequencies
        error('ACC_Unregularized:DimensionMismatch', ...
            'Bright- and dark-zone ATFs must share loudspeaker/frequency dimensions.');
    end

    nearly_singular_state = warning('query', ...
        'MATLAB:nearlySingularMatrix');
    singular_state = warning('query', 'MATLAB:singularMatrix');
    warning_cleanup = onCleanup(@() local_restore_warnings( ...
        nearly_singular_state, singular_state));
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');

    w = zeros(loudspeakers, num_frequencies, 'like', 1i);
    for frequency_index = 1:num_frequencies
        HB = HB_ctrl(:, :, frequency_index);
        HD = HD_ctrl(:, :, frequency_index);
        RB = HB' * HB;
        RD = HD' * HD;
        quotient = RD \ RB;
        [~, w(:, frequency_index)] = MaxEigenvector(quotient);
    end
end

function local_restore_warnings(nearly_singular_state, singular_state)
    warning(nearly_singular_state);
    warning(singular_state);
end
