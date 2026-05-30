function AE = calculate_AE(w, H_B, virtual_source_idx)
% CALCULATE_AE - Computes the Array Effort (AE) in dB.
    
    w0 = zeros(size(w));
    w0(virtual_source_idx) = 1; % reference
    
    % Add a small epsilon for numerical stability
    epsilon = 1e-12;
    
    effort_ratio = w' * w;
    RB = (H_B' * H_B);
    energy_ratio = (w0' * RB * w0) / (w' * RB * w + epsilon);
% 单独发声的虚拟声源声压作为基准值，然后计算出这个阵列在当前控制w下
% 在该声场需要产生的同等声压的增益比，然后来计算出w能量应该放大多少倍，
% 然后放大w，此时的阵列增益叫做AE
    
    AE = 10 * log10(effort_ratio * energy_ratio + epsilon);
end