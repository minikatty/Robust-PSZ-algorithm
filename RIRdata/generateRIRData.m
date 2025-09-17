%% array
clc;close;
clear;
addpath(genpath(pwd));
load("array.mat");
layout_show(array);

%% generate data 
% setup
beta = 0.3;
frebin = 20; % spectral resolution, frequency bin
f = frebin:frebin:4000;
c = 343; % about 17.65 degrees Celsius, ((342/331.45)^2 - 1)*273
k = 2*pi*f/c; % wave number
n = 2*length(f);
fs = frebin*n;

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
for i = 1:size(s, 1)
    parfor j = 1:size(bCtrPtsPositions, 1)
        rir = rir_generator(c, fs, bCtrPtsPositions(j, :), s(i, :), roomSize, beta, n);
        tempH = fft(rir, n);
        HB(j, i, :) = 2*tempH(2:floor(n/2)+1); % Discard dc component
        
        rir = rir_generator(c, fs, bPerPtsPositions(j, :), s(i, :), roomSize, beta, n);
        tempH = fft(rir, n);
        HBE(j, i, :) = 2*tempH(2:floor(n/2)+1); 

        rir = rir_generator(c, fs, dCtrPtsPositions(j, :), s(i, :), roomSize, beta, n);
        tempH = fft(rir, n);
        HD(j, i, :) = 2*tempH(2:floor(n/2)+1); 

        rir = rir_generator(c, fs, dPerPtsPositions(j, :), s(i, :), roomSize, beta, n);
        tempH = fft(rir, n);
        HDE(j, i, :) = 2*tempH(2:floor(n/2)+1); 
    end
end
% generate RIR adding gaussian noisy 

snr = 15:(25-15)/(NumberNAdd-1):25; 
% generate different snr data, space regard the number of adding noise

for ij = 1:NumberNAdd
    for i = 1:size(s, 1)
        parfor j = 1:size(bCtrPtsPositions, 1)
            rir0 = rir_generator(c, fs, bCtrPtsPositions(j, :), s(i, :), roomSize, beta, n);
            rir = awgn(rir0, snr(ij), 'measured');
            tempH = fft(rir, n);
            HBMeasured(j, i, :, ij) = 2*tempH(2:floor(n/2)+1);
    
            rir0 = rir_generator(c, fs, dCtrPtsPositions(j, :), s(i, :), roomSize, beta, n);
            rir = awgn(rir0, snr(ij), 'measured');
            tempH = fft(rir, n);
            HDMeasured(j, i, :, ij) = 2*tempH(2:floor(n/2)+1); 
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
save RIRdata/HDMeasured.mat HDMeasured;
save RIRdata/para.mat para;