function[w] = ACC_PM(HB_ctrl, HD_ctrl, HB_Desired, para) 
% J = kappa * P_d^H * P_d ...
%     + (1 - kappa) * (P_b - P_b_hat)^H * (P_b - P_b_hat);
    kappa = para.kappa;
    target_freqs = para.target_freqs;
    num_frePoint = length(target_freqs);
    L  = size(HB_ctrl, 2);
    w = zeros(L, num_frePoint);
%     freq = para.freq;
    for i = 1:num_frePoint
        HB = squeeze(HB_ctrl(:, :, i));
        HD = squeeze(HD_ctrl(:, :, i));
        pBd = HB_Desired(:, i);
        R = kappa*(HD'*HD) + (1-kappa)*(HB'*HB);
        belta = eps;%max(eig(R))*(1e-5);
        % belta = 1e-5*max(eig(R));%max(eig(R))*(2e-5);
        w(:, i) = (R + belta*eye(size(HB, 2))) \ ((1-kappa) * (HB' * pBd));  
    end
end
