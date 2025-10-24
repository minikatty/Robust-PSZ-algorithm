function freq_params = configure_freq_parameters(varargin)
% CONFIGURE_FREQ_PARAMETERS - Configures and computes ESSENTIAL frequency-domain parameters.
%
% This function centralizes the definition of sampling rates, FFT settings,
% and target frequency mappings for the entire project, ensuring consistency.
%
% It returns a lean structure containing only the parameters necessary for
% downstream analysis.
%
% Outputs:
%   freq_params - A structure containing all relevant parameters.
%
% Example:
%   freq_params = configure_freq_parameters();
%   freq_params = configure_freq_parameters('fs_downsampled', 8000, 'target_freq_end', 5000);

    p = inputParser;
    addParameter(p, 'fs_native', 48000, @isnumeric);
    addParameter(p, 'fs_downsampled', 16000, @isnumeric);
    addParameter(p, 'rir_duration_ms', 128, @isnumeric);
    addParameter(p, 'target_freq_start', 100, @isnumeric);
    addParameter(p, 'target_freq_end', 4000, @isnumeric);
    addParameter(p, 'target_freq_step', 25, @isnumeric);
    parse(p, varargin{:});

    params = p.Results;

    % --- Core Calculations ---
    rir_len_native = floor(params.rir_duration_ms / 1e3 * params.fs_native);
    [p_resample, q_resample] = rat(params.fs_downsampled / params.fs_native);
    nfft = round(rir_len_native * p_resample / q_resample);
    
    num_freq_bins = floor(nfft / 2) + 1;
    freq_axis = (0:num_freq_bins-1)' * params.fs_downsampled / nfft;
    fft_freq_resolution = params.fs_downsampled / nfft;

    f_target_physical = (params.target_freq_start : params.target_freq_step : params.target_freq_end)';
    f_target_indices = round(f_target_physical / fft_freq_resolution) + 1;
    
    % Ensure indices are within bounds
    f_target_indices(f_target_indices > num_freq_bins) = num_freq_bins;
    
    f_target_actual = freq_axis(f_target_indices);

    % ---the essential parameters into the output structure ---
    freq_params.fs = params.fs_downsampled; 
    freq_params.nfft = nfft;
    freq_params.freq_axis = freq_axis;
    freq_params.target_indices = f_target_indices; % for plot
    freq_params.target_freqs = f_target_actual; 
    freq_params.p_resample = p_resample;
    freq_params.q_resample = q_resample;
    freq_params.rir_len_native = rir_len_native;
    freq_params.num_freq_bins = num_freq_bins;
    
    % Optional: Display summary when called directly for debugging
    if nargout == 0
        fprintf('--- Frequency Parameters Summary ---\n');
        disp(freq_params);
    end
end

