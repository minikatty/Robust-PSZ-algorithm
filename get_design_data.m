clc; clear; close all;

addpath(genpath('data/'));
current_mode = 'temperature';
data_dir = fullfile('data', 'SimulateRIR', current_mode);
numFiles = length(20.00:0.10:25);
for design_env_para = 20.00:0.10:25
    design_filename = get_data_filename(data_dir, current_mode, design_env_para);
    data = load(design_filename);
    current_BZ = data.ATF_BZ.ctrl;
    current_DZ = data.ATF_DZ.ctrl;
    if design_env_para == 20.00
        sum_ATF_BZ = zeros(size(current_BZ));
        sum_ATF_DZ = zeros(size(current_DZ));
    end
    sum_ATF_BZ = sum_ATF_BZ + current_BZ;
    sum_ATF_DZ = sum_ATF_DZ + current_DZ;
end
    
Mean_Data.BZ = sum_ATF_BZ / numFiles;
Mean_Data.DZ = sum_ATF_DZ / numFiles;


figure;
subplot(3,1,1);
% 假设数据是频域数据，画幅值
plot(20*log10(abs(squeeze(Mean_Data.BZ(1,1,:))))); 
title('Mean ATF BZ (Magnitude)');
ylabel('Magnitude (dB)'); grid on;

subplot(3,1,2);
plot(20*log10(abs(squeeze(Mean_Data.BZ(1,1,:)))));
title('Mean ATF DZ (Magnitude)');
ylabel('Magnitude (dB)'); grid on;

subplot(3,1,3);
plot(20*log10(abs(squeeze(current_BZ(1,1,:)))));
title('Mean Current BZ(Magnitude)');
ylabel('Magnitude (dB)'); grid on;
% 5. 保存结果
save('mean_design_data.mat', 'Mean_Data');