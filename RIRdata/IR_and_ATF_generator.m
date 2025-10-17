%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%       Generate RIR and ATFs data 
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% array
clc;
% close;
clear;
addpath(genpath(pwd));
load("RIRdata/array.mat");
layout_show(array); % plot the layout and scenes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters settings
beta = 0.3;           % T60 reverberation time, all walls are same
% target freqz band [100,4e3]<====>idx:[5:end]
c = 343;              % sound speed (m/s)
fs = 48000;           % sample rate (Hz)
% nfft = fs/frebin;
% nfft = 2*length(f) + 1; 
% odd number get a symmetric spectrum, no Nyquist freqz
% nf = floor(nfft/2) + 1; % single sided freqz number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% truncated length and NFFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nsamples = ceil(T60 * fs); % nsamples
% nsamples = nfft;
truncatad_time = 128; % unit: ms, include early reverberation
len_truncated = floor(truncatad_time/1e3 * fs); 
% IR truncated length: 50ms include early reflection
rir_len = len_truncated;  % two-side spectrum
f_ds = 16000;
[p, q] = rat(f_ds/fs);
frebin = 1e3/truncatad_time ; % frequency resolution for plot
target_f_start = 200;   % target control freq band  100~4000   
target_f_end = 4000;
frebin_tar = 25;
f = target_f_start:frebin_tar:target_f_end; % desired control freq resolution
target_f = ceil(f./frebin); % target f idx
n_tarf = length(target_f); % number of positive freq points, no dc component
% rir_len_eval = rir_len; % for evalute the input speech/music
nf_sing = floor(rir_len*f_ds/(2*fs))+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              array layout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s  = array.s; %
bCtrPtsPositions = array.bCtrPtsPositions; 
dCtrPtsPositions = array.dCtrPtsPositions;
bPerPtsPositions = array.bPerPtsPositions;
dPerPtsPositions = array.dPerPtsPositions;
roomSize = array.roomSize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SNR  range
NumberNAdd = 51;        % Number of noise additions
snr_range = [20 40];   %  30-50 is ok, can try different
snr = linspace(snr_range(1), snr_range(end), NumberNAdd);

% tmeperature


% position mismatch
% only 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  initialization of matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nCtr = size(bCtrPtsPositions, 1); % 控制点个数: 96(48 control, 48 monitor)
                                  % same as evaluate points
nSrc = size(s, 1);                % 源个数: 48
% nCtr = 2;nSrc = 2;                % debug
HBMeasured = zeros(nCtr, nSrc, nf_sing, NumberNAdd);
HDMeasured = zeros(nCtr, nSrc, nf_sing, NumberNAdd);
% [mic,lsk,n_freq,mearsurec_counts]
HB_eval = HBMeasured;
HD_eval = HDMeasured;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Noisy RIR generator
parpool('local', 8);
rng(2025); % for reproductivity
for ij = 1:NumberNAdd
    HBtemp = zeros(nCtr, nSrc, nf_sing); % ATFs temp memery
    HDtemp = zeros(nCtr, nSrc, nf_sing);
    HBtemp_eval = HBtemp;
    HDtemp_eval = HDtemp;
    IR_BZ_ctrl_tmp = zeros(nCtr, nSrc, round(rir_len*p/q));
    IR_BZ_eval_tmp = IR_BZ_ctrl_tmp;
    IR_DZ_ctrl_tmp = IR_BZ_ctrl_tmp;
    IR_DZ_eval_tmp = IR_BZ_ctrl_tmp;
    snr_ij = snr(ij);
   for i = 1:nSrc
       src_pos_i = s(i,:);
        parfor j = 1:nCtr % parfor
            %====================BZ=========================================
            rir_ctrl = rir_generator(c, fs, bCtrPtsPositions(j,:), src_pos_i, ...
                roomSize, beta, rir_len);
            rir_ctrl = resample(rir_ctrl, p, q);
            rir_eval = rir_generator(c, fs, bPerPtsPositions(j,:), src_pos_i, ...
                roomSize, beta, rir_len);
            rir_eval = resample(rir_eval, p, q); % downsample
            rir_ctrl = awgn(rir_ctrl, snr_ij, 'measured'); % noisy rirs
            rir_eval = awgn(rir_eval, snr_ij, 'measured');
            IR_BZ_ctrl_tmp(j,i,:) = rir_ctrl;
            IR_BZ_eval_tmp(j,i,:) = rir_eval;
            tempH = fft(rir_ctrl);
            HBtemp(j, i, :) = tempH(1:nf_sing); % single side spectrum
            tempH = fft(rir_eval); % reuse
            HBtemp_eval(j, i, :) = tempH(1:nf_sing);
            %==========================DZ==================================
            rir_ctrl = rir_generator(c, fs, dCtrPtsPositions(j,:), src_pos_i, ...
                roomSize, beta, rir_len);
            rir_ctrl = resample(rir_ctrl, p, q);
            rir_eval = rir_generator(c, fs, dPerPtsPositions(j,:), src_pos_i, ...
                roomSize, beta, rir_len);
            rir_eval = resample(rir_eval, p, q);
            rir_ctrl = awgn(rir_ctrl, snr_ij, 'measured');
            rir_eval = awgn(rir_eval, snr_ij, 'measured');
            IR_DZ_ctrl_tmp(j,i,:) = rir_ctrl;
            IR_DZ_eval_tmp(j,i,:) = rir_eval;
            tempH = fft(rir_ctrl);
            HDtemp(j, i, :) = tempH(1:nf_sing);
            tempH = fft(rir_eval);
            HDtemp_eval(j, i, :) = tempH(1:nf_sing);
        end
   end
    IR.IR_BZ_ctrl(:, :, :, ij) = IR_BZ_ctrl_tmp;
    IR.IR_BZ_eval(:, :, :, ij) = IR_BZ_eval_tmp;
    IR.IR_DZ_ctrl(:, :, :, ij) = IR_DZ_ctrl_tmp;
    IR.IR_DZ_eval(:, :, :, ij) = IR_DZ_eval_tmp;
    HBMeasured(:, :, :, ij) = HBtemp;
    HDMeasured(:, :, :, ij) = HDtemp;
    HB_eval(:, :, :, ij) = HBtemp_eval;
    HD_eval(:, :, :, ij) = HDtemp_eval;
end
% nfft = f_ds*truncatad_time/1e3;
delete(gcp('nocreate'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save data
para.f = (0:nf_sing-1) * frebin;   % all freq idx of single sided spectrum
para.ftar = para.f(target_f+1); % true control freq
para.fs = f_ds;
para.dimHBMeasured = '[nCtr, nSrc, nFreq, nNoise]';
% IR
para.snr = snr;

% save('RIRdata/HB1004.mat', '-v7.3', 'HBMeasured', 'HB_eval');
% save('RIRdata/HD1004.mat', '-v7.3', 'HDMeasured', 'HD_eval');
save('RIRdata/IR_data1004.mat', '-v7.3','IR');
save('RIRdata/para1004.mat', 'para');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%