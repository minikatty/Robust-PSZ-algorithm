clear;clc;close all;
algorithms_to_test = {'RACC_PM_GLS','RACC_PM_Subpro'};
current_mode = 'temperature';
data_dir = fullfile('data', 'SimulateRIR', current_mode);
para_file = fullfile(data_dir, 'para.mat');

load(para_file); % Loads 'para' setting file

freq_params = para.freq_params; 
target_freqs = freq_params.target_freqs;
num_target_freqs = length(target_freqs);
design_env_para = 22.5; % temperature
design_filename = get_data_filename(data_dir,current_mode,design_env_para);
design_data = load(design_filename);

num_points = 15;
rho_list = logspace(-3, 3, num_points); 
plane_wav_flag = 1;
if plane_wav_flag
    load('./data/ATF_desired_plane.mat'); 
    ATF_desired = ATF_desired_plane;
else
    ATF_desired = squeeze(ATF_BZ_ctrl(:, virtual_src_idx, :)); 
    % use nominal env
end
design_data.ATF_desired = ATF_desired;

[gamma, epsilon] = get_bound_paras(current_mode, design_filename);
freq_params.epsilon = epsilon;
freq_params.gamma = gamma;
freq_params.scale = 1e-3;
freq_params.virtual_src_idx = 13; % target sound field,PM-related
% ACC_PM:
kappa = 0.5; % the weight for BZ
freq_params.kappa = kappa;
% POTDC_RACC/wcRACC/RPM: bound paras above
% RACC_PM_GLS/RACC_PM_Subpro:
% alpha for RACC-PM 

ATF_BZ_ctrl = design_data.ATF_BZ.ctrl;ATF_DZ_ctrl = design_data.ATF_DZ.ctrl;
ATF_BZ_eval = design_data.ATF_BZ.eval;ATF_DZ_eval = design_data.ATF_DZ.eval; %eval,ctrl
wopt = wcACC(ATF_BZ_ctrl, ATF_DZ_ctrl,freq_params);
if isfield(wopt,'scale')
    scale = wopt.scale;
    filters_w = wopt.w;
else
    scale = nan;
    filters_w = wopt;
end
margin = 2; % Leaving a 3 dB margin to the inaccuracies upper bound from wcRACC.
alpha_AC = zeros(num_target_freqs,1);

for i = 1:num_target_freqs
% Extract frequency-dependent data for this iteration.
    w_f   = filters_w(:, i);
    H_B_f = ATF_BZ_eval(:, :, i); 
    H_D_f = ATF_DZ_eval(:, :, i);
    alpha_AC(i) = 10^((real(calculate_AC(w_f, H_B_f, H_D_f))-margin)/10);
end
freq_params.mu = 1; % equal for RPM-ACC & PM
freq_params.alpha = alpha_AC; % refer to wcRACC-2dB

obj_freq = [500, 1000, 4000]; % 
num_freqs = length(obj_freq);
freq_idx = [9, 19, 79]; % corresonding to obj_freq

w_opt_RACC = design_filters('wcRACC',design_data,freq_params);
w_opt_RPM = design_filters('RPM',design_data,freq_params);

% Allocate Memery [num_freqs x num_points]
AC_Sub_res   = zeros(num_freqs, num_points);
NSRE_Sub_res = zeros(num_freqs, num_points);
AC_GLS_res   = zeros(num_freqs, num_points);
NSRE_GLS_res = zeros(num_freqs, num_points);

% base_method metric
AC_RACC_val  = zeros(num_freqs, 1);
AC_RPM_val   = zeros(num_freqs, 1);
NSRE_RPM_val = zeros(num_freqs, 1);

for f_idx = 1:num_freqs
    
    c_freq = obj_freq(f_idx);
    f_id   = freq_idx(f_idx); 
    
    H_B_f = squeeze(ATF_BZ_eval(:, :, f_id));
    H_D_f = squeeze(ATF_DZ_eval(:, :, f_id));
    p_d_f = ATF_desired(:, f_id); 
    NSRE_H_B_f =  squeeze(ATF_BZ_ctrl(:, :, f_id));
    
    w_RACC_f = w_opt_RACC(:, f_id);
    w_RPM_f  = w_opt_RPM(:, f_id);
    
    AC_RACC_val(f_idx)  = real(calculate_AC(w_RACC_f, H_B_f, H_D_f));
    AC_RPM_val(f_idx)   = real(calculate_AC(w_RPM_f,  H_B_f, H_D_f));
    NSRE_RPM_val(f_idx) = calculate_NSRE(w_RPM_f, NSRE_H_B_f, p_d_f);
 
  

    
    fprintf('\n======== 开始处理频点: %d Hz ========\n', c_freq);

    for k = 1:num_points
        current_rho = rho_list(k);
        fprintf('  >> 计算 rho = %.2e (%d/%d)...\n', current_rho, k, num_points);
        
        freq_params.rho = current_rho; 
        freq_params.para_gamma = freq_params.mu + current_rho * freq_params.alpha;
        
        wopt_Sub(f_idx,k)= RACC_PM_Sub(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, freq_params, f_id);
        
        wopt_GLS(f_idx,k) = RACC_PM_GLS(ATF_BZ_ctrl, ATF_DZ_ctrl, ATF_desired, freq_params, f_id);
        
        w_Sub_f = wopt_Sub(f_idx,k).w(:, f_id);
        w_GLS_f = wopt_GLS(f_idx,k).w(:, f_id);
        
        AC_Sub_res(f_idx, k)   = real(calculate_AC(w_Sub_f, H_B_f, H_D_f));
        NSRE_Sub_res(f_idx, k) = calculate_NSRE(w_Sub_f, NSRE_H_B_f, p_d_f);
        
        AC_GLS_res(f_idx, k)   = real(calculate_AC(w_GLS_f, H_B_f, H_D_f));
        NSRE_GLS_res(f_idx, k) = calculate_NSRE(w_GLS_f, NSRE_H_B_f, p_d_f);
    end      
end
timestamp_str = datetime('now', 'format','yyyyMMdd');    
% joint filename
save_filename = "Pareto_Tradeoff_Results346_" + string(timestamp_str) + ".mat";
save(save_filename, ...
    ... % basic settings 
    'rho_list', 'num_points', ... 
    'obj_freq', 'freq_idx', 'num_freqs', ...
    ...
    ... % metric matrix 'AC_Sub_res', 'NSRE_Sub_res', ...
    'AC_GLS_res', 'NSRE_GLS_res', ...
    ...
    'w_GLS_f');
    % ... % benchmark metric
     % 'NSRE_RPM_val', ... %'AC_RACC_val', 'AC_RPM_val', ...... % filters
      % 'w_Sub_f','w_opt_RPM', 'w_opt_RACC',
fprintf('数据保存成功！\n');
%% plot
% color
c_RACC = [0.4 0.7 0.4]; 
c_RPM  = [0.9 0.6 0.1]; 
c_GLS  = [0.8 0.3 0.3]; 
c_Sub  = [0.2 0.6 0.8]; 

for f_idx = 1:num_freqs
    current_freq = obj_freq(f_idx);

    % extract current idx
    ac_rpm   = AC_RPM_val(f_idx);
    ac_racc  = AC_RACC_val(f_idx);
    nre_rpm  = NSRE_RPM_val(f_idx);

    ac_gls   = AC_GLS_res(f_idx, :);
    ac_sub   = AC_Sub_res(f_idx, :);
    nre_gls  = NSRE_GLS_res(f_idx, :);
    nre_sub  = NSRE_Sub_res(f_idx, :);
    % =========================================================================
    % AC Plot
    % =========================================================================
    fig_AC = figure('Position',[100, 200, 520, 420], 'Color', 'w', 'Name', sprintf('AC_%dHz', current_freq));
    hold on; grid on; box on;

    h_RACC_AC = plot([rho_list(1), rho_list(end)],[ac_racc, ac_racc], '--', 'Color', c_RACC, 'LineWidth', 2.5);
    h_RPM_AC  = plot([rho_list(1), rho_list(end)],[ac_rpm, ac_rpm], '-.', 'Color', c_RPM, 'LineWidth', 2.5);

    h_GLS_AC = plot(rho_list, ac_gls, '-s', 'Color', c_GLS, 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'w');
    h_Sub_AC = plot(rho_list, ac_sub, '-o', 'Color', c_Sub, 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'w');

    set(gca, 'XScale', 'log', 'FontSize', 12, 'FontName', 'Times New Roman', 'LineWidth', 1.2);
    xlabel('Weighting Parameter \rho', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('AC: [dB]', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([rho_list(1), rho_list(end)]);
    title(sprintf('Acoustic Contrast at %d Hz', current_freq), 'FontSize', 14);

    legend([h_RACC_AC, h_GLS_AC, h_Sub_AC, h_RPM_AC], ...
        {'wcRACC', 'Proposed wcRACC-PM', 'Proposed wcRACC-PM-sub', 'Proposed RPM'}, ...
        'Location', 'best', 'FontSize', 11, 'FontName', 'Times New Roman');
    % =========================================================================
    %  NRE Plot
    % =========================================================================
    fig_NRE = figure('Position',[650, 200, 520, 420], 'Color', 'w', 'Name', sprintf('NRE_%dHz', current_freq));
    hold on; grid on; box on;

    h_RPM_NRE  = plot([rho_list(1), rho_list(end)], [nre_rpm, nre_rpm], '--', 'Color', c_RPM, 'LineWidth', 2.5);

    h_GLS_NRE = plot(rho_list, nre_gls, '-s', 'Color', c_GLS, 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'w');
    h_Sub_NRE = plot(rho_list, nre_sub, '-o', 'Color', c_Sub, 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'w');

    set(gca, 'XScale', 'log', 'FontSize', 12, 'FontName', 'Times New Roman', 'LineWidth', 1.2);
    xlabel('Weighting Parameter \rho', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('NRE: [dB]', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([rho_list(1), rho_list(end)]);
    title(sprintf('Reproduction Error at %d Hz', current_freq), 'FontSize', 14);

    legend([h_GLS_NRE, h_Sub_NRE, h_RPM_NRE], ...
        {'Proposed wcRACC-PM', 'Proposed RACC-PM-sub', 'Proposed RPM'}, ...
        'Location', 'best', 'FontSize', 11, 'FontName', 'Times New Roman');
end