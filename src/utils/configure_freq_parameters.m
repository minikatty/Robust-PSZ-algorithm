function freq_params = configure_freq_parameters(varargin)
% CONFIGURE_FREQ_PARAMETERS - Configures and computes ESSENTIAL frequency-domain parameters.
%
% This function centralizes the definition of sampling rates, FFT settings,
% and target frequency mappings for the entire project, ensuring consistency.
%
% It returns a lean structure containing only the parameters necessary for analysis.
%
% Outputs:
%   freq_params - A structure containing all relevant parameters.
%
% Example:
%   freq_params = configure_freq_parameters();
%   freq_params = configure_freq_parameters('target_freq_end', 5000);

p = inputParser;
addParameter(p, 'fs', 16000, @isnumeric);
addParameter(p, 'rir_duration', 150, @isnumeric); % ms
addParameter(p, 'target_freq_start', 100, @isnumeric);
addParameter(p, 'target_freq_end', 8000, @isnumeric);
addParameter(p, 'target_freq_step', 40, @isnumeric);
parse(p, varargin{:});
freq_params = p.Results;
rir_duration = freq_params.rir_duration;
fs = freq_params.fs;
target_freq_step = freq_params.target_freq_step;
target_freq_start = freq_params.target_freq_start;
target_freq_end = freq_params.target_freq_end;

% --- Core Calculations ---
rir_len = floor(rir_duration/ 1e3 * fs);
nfft = rir_len;
fft_freq_resolution = fs / nfft;

num_freq_bins = floor(nfft / 2) + 1;
freq_axis = (0:num_freq_bins-1)'*fft_freq_resolution;

f_target = unique([target_freq_start : target_freq_step : target_freq_end, target_freq_end])';
f_target_indices = round(f_target / fft_freq_resolution) + 1;

% Ensure indices are within bounds
f_target_indices(f_target_indices > num_freq_bins) = num_freq_bins;

f_target_actual = freq_axis(f_target_indices);

% ---the essential parameters into the output structure ---
% freq_params.freq_axis = freq_axis;
% freq_params.target_indices = f_target_indices; % for plot
freq_params.target_freqs = f_target_actual; 
freq_params.rir_len = rir_len;
freq_params.fs = fs;
freq_params.target_freq_step = target_freq_step;
freq_params.target_freq_start = target_freq_start;
freq_params.target_freq_end = target_freq_end;
freq_params.nfft = nfft;
freq_params.f_target_indices = f_target_indices;

end
