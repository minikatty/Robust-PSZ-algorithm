%% array
clc;close;
clear;
addpath(genpath(pwd));
load("array.mat");
% layout_show(array); % plot the layout and scenes

%% generate data 
% setup
beta = 0.3; %
frebin = 20; % spectral resolution, frequency bin
f = frebin:frebin:4000; % freq_bands 
c = 343; % about 17.65 degrees Celsius, ((342/331.45)^2 - 1)*273
k = 2*pi*f/c; % wave number
fs = 16000;% sample frq 
T60 = beta;
nsamples = ceil(T60 * fs); % nsamples
nfft = 2*length(f);
len_truncated = ceil(50/1000 * fs); % 50 ms----ITU-T G.167, G.168 suggest,reserve the early reviberation
% for speech enhancement

s = array.s;
bCtrPtsPositions = array.bCtrPtsPositions;
bPerPtsPositions = array.bPerPtsPositions; % 明区评测点
dCtrPtsPositions = array.dCtrPtsPositions;
dPerPtsPositions = array.dPerPtsPositions; % 暗区评测点

roomSize = array.roomSize;
HB = zeros(size(bCtrPtsPositions, 1), size(s, 1), length(f));
HD = HB;
HBE = HB;
HDE = HB;
NumberNAdd = 1;% Number of noise additions
HBMeasured = zeros(size(HB, 1), size(HB, 2), size(HB, 3), NumberNAdd);
HDMeasured = HBMeasured;
% generate true and evaluation RIR using toolbox
% for i = 1:size(s, 1)
%     parfor j = 1:size(bCtrPtsPositions, 1)
%         rir = rir_generator(c, fs, bCtrPtsPositions(j, :), s(i, :), roomSize, beta, n);
%         tempH = fft(rir, n);
%         HB(j, i, :) = 2*tempH(2:floor(n/2)+1); % Discard dc component
% 
%         rir = rir_generator(c, fs, bPerPtsPositions(j, :), s(i, :), roomSize, beta, n);
%         tempH = fft(rir, n);
%         HBE(j, i, :) = 2*tempH(2:floor(n/2)+1); 
% 
%         rir = rir_generator(c, fs, dCtrPtsPositions(j, :), s(i, :), roomSize, beta, n);
%         tempH = fft(rir, n);
%         HD(j, i, :) = 2*tempH(2:floor(n/2)+1); 
% 
%         rir = rir_generator(c, fs, dPerPtsPositions(j, :), s(i, :), roomSize, beta, n);
%         tempH = fft(rir, n);
%         HDE(j, i, :) = 2*tempH(2:floor(n/2)+1); 
%     end
% end
% generate RIR adding gaussian noisy 
snr_range = [30 50];
snr = snr_range(1):(snr_range(end)-snr_range(1))/(NumberNAdd-1)...
    :snr_range(end); 
% generate different snr data, space regard the number of adding noise

for ij = 1:NumberNAdd
    for i = 1:size(s, 1)
        for j = 1:size(bCtrPtsPositions, 1) % parfor improve efficiency remove for debug
            rir0 = rir_generator(c, fs, bCtrPtsPositions(j, :), s(i, :), roomSize, beta, nsamples);
            % delay_sec = 1.2*norm(bCtrPtsPositions(j, :) - s(i, :)) / c; %
            % direct time,not goog enough for pratical auditory experience
            % len_truncated = round(delay_sec * fs);
            plot_rir(rir0,fs)
            rir_truncate = rir0(1:len_truncated); % truncate the rir for room compensation
            plot_rir(rir_truncate,fs)
            rir = awgn(rir_truncate, snr(ij), 'measured');
            tempH = fft(rir, nfft);
            HBMeasured(j, i, :, ij) = 2*tempH(2:floor(nfft/2)+1);
            % rir0 is tmp variable
            rir0 = rir_generator(c, fs, dCtrPtsPositions(j, :), s(i, :), roomSize, beta, nsamples);
            rir = awgn(rir0, snr(ij), 'measured');
            tempH = fft(rir, nfft);
            HDMeasured(j, i, :, ij) = 2*tempH(2:floor(nfft/2)+1); 
        end
    end
end
para.f = f;
% save RIR
delete(gcp('nocreate'));
save('RIRdata/HB.mat', 'HB');
save RIRdata/HBE.mat HBE;
save RIRdata/HD.mat HD;
save RIRdata/HDE.mat HDE;
save RIRdata/HBMeasured.mat HBMeasured; 
%[number_points, number_lsk, number_freq, number_add]
save RIRdata/HDMeasured.mat HDMeasured;
save RIRdata/para.mat para;