function [AE, w] = calculate_AE(w, H_B, virtual_source_idx, target_SPL)
% compute the corresponding Array Effort after calibrating the SPL
%
% imput:
%   w      : 算法算出的原始权重向量 (N x 1)
%   H_B   : 亮区麦克风的传递函数矩阵 (M_bright x N)
%   target_SPL : 目标声压级，例如 76 (dB SPL)
%
% 输出:
%   AE: 标定后的阵列增益 (dB)
%   w : 标定后的物理权重 
    w0 = zeros(size(w));
    w0(virtual_source_idx) = 1; % reference
    RB = H_B' * H_B;
    P_ref = 20e-6;  % ref pressure in the air
    P_target_Pa = P_ref * db2mag(target_SPL); 
    p_raw = H_B * w; % M*1
    p_current_RMS = sqrt(mean(abs(p_raw).^2)) + eps; % pressure density
    alpha = P_target_Pa / p_current_RMS; % scaling factor
    w = w * alpha; % calibrate the filter weights amplitude    
    % 单独发声的虚拟声源声压作为基准值，然后计算出这个阵列在当前控制w下
    % 在该声场需要产生的同等声压的增益比，然后来计算出w能量应该放大多少倍，
    % 然后放大w，此时的阵列增益叫做AE   
    energy_ratio = real((w0' * RB * w0) / (w' * RB * w + eps)); 
    % energy_ratio is 1/E_{ref} in the paper
    effort_ratio = w' * w;
    AE = 10 * log10(effort_ratio * energy_ratio + eps);
end


