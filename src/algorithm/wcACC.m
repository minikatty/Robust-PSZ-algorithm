function res = wcACC(HB_ctrl, HD_ctrl, wcACC_para)
    frePoint = wcACC_para.target_freqs;
    % gammaB = wcACC_para.gamma.B;
    % gammaD = wcACC_para.gamma.D;
    num_frePoint = length(frePoint);
    L  = size(HB_ctrl, 2);
    w = zeros(L,num_frePoint);
    scale1 = wcACC_para.scale; % 1
    % scale2 = wcACC_para.scale; % for simulation:1e-2, real measured: 1e-4
    for i = 1:num_frePoint
        HB = squeeze(HB_ctrl(:, :, i));
        HD = squeeze(HD_ctrl(:, :, i));
        RB = HB'*HB;
        RD = HD'*HD;        
        scale = scale1;
        GammaB = 0.75 * norm(RB,'fro') /100; % refer to paper
        GammaD = GammaB;        
        A = (RD + GammaD *eye(size(RD))) \ (RB - GammaB *eye(size(RB)));
        [~, w(:, i)] = MaxEigenvector(A);        
    end
    res.w =w;
    res.scale = scale;
end

