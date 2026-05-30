% %% data processing
clc; clear; close all;
data_dir = fullfile('data', 'Cabin_Measurements');
addpath(data_dir);
addpath('src');
% 
init_file = fullfile(data_dir, 'ir_zhenren_1.mat');
ir_init = load(init_file); 
fs = ir_init.fs_ir;
lsk_num = size(ir_init.saved_data, 1);
[~, mic_num] = size(ir_init.saved_data{1});

freq_params = configure_freq_parameters('fs', fs, 'rir_duration', 100, ...
    'target_freq_start', 100, 'target_freq_end', 5000,'target_freq_step', 100); 

ir_length = freq_params.rir_len;
BZ_idx = 1:4; 
DZ_idx = 5:8;

% --- Preallocation)  ---
ir_mtx = zeros(mic_num, lsk_num, ir_length); 
Nominal_ir = zeros(mic_num, lsk_num, ir_length); 

%  Fade-out window
fade_time = 10; % ms
fade_samples = round(fade_time / 1000 * fs); 
w_tail = hann(2 * fade_samples);
w_tail = w_tail(fade_samples+1 : end);

w_tail_3d = reshape(w_tail, 1, 1, fade_samples); 


file_num = 60;
for k = 1 : file_num
    filename = fullfile(data_dir, sprintf('ir_zhenren_%d.mat', k));
    ir = load(filename); 
    saved_data = ir.saved_data;

    for i = 1:lsk_num
        ir_mtx(:, i, :) = saved_data{i}(1:ir_length, :)';
    end 

    ir_mtx(:, :, end-fade_samples+1:end) = ir_mtx(:, :, end-fade_samples+1:end) .* w_tail_3d;

    Nominal_ir = Nominal_ir + ir_mtx;

    ATF_tmp = compute_atf(ir_mtx, freq_params);
    ATF.ATF_BZ = ATF_tmp(BZ_idx, :, :);
    ATF.ATF_DZ = ATF_tmp(DZ_idx, :, :);    

    saved_filename = fullfile(data_dir, sprintf('ATF_%d.mat', k));
    save(saved_filename, '-struct', 'ATF');    
end

Nominal_ATF.Nominal_ir = Nominal_ir ./ file_num;
ATF_tmp = compute_atf(Nominal_ir, freq_params);

Nominal_ATF.ATF_BZ = ATF_tmp(BZ_idx, :, :);
Nominal_ATF.ATF_DZ = ATF_tmp(DZ_idx, :, :);
Nominal_ATF.freq_params = freq_params;

nominal_saved_filename = fullfile(data_dir, 'Nominal_ATF.mat');
save(nominal_saved_filename, '-struct', 'Nominal_ATF');

disp('=== 所有的 ATF 与 Nominal ATF 数据已处理并保存完毕！ ===');