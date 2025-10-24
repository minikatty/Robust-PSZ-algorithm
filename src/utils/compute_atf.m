function atf_single_sided = compute_atf(rir_matrix, freq_params)
% COMPUTE_ATF - Converts time-domain RIR to frequency-domain ATF (single-sided).
%
% Usage:
%   atf = compute_atf(rir_matrix, freq_params);
%
% Inputs:
%   rir_matrix  - Time-domain signal [M x N x L]
%   freq_params - Structure containing frequency parameters (must have .nfft)
%
% Outputs:
%   atf_single_sided - Single-sided frequency spectrum [M x N x (nfft/2+1)]

    % 1. 获取 FFT 参数
    nfft = freq_params.nfft;
    
    % 2. 计算单边谱所需的频点数
    % (freq_params 可能有这个字段，但重新计算一次也很简单且安全)
    num_freq_bins = freq_params.num_freq_bins;
    
    % 3. 执行 FFT (沿第3维，即时间轴)
    atf_full = fft(rir_matrix, nfft, 3);
    
    % 4. 截取单边谱
    atf_single_sided = atf_full(:, :, 1:num_freq_bins);

end