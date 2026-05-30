function[w] = PM(HB_ctrl, HD_ctrl, HB_Desired)
    % PM algorithm implementation  
    [~, L, num_frePoint]  = size(HB_ctrl);
    % target_pressure =
    w = zeros(L,num_frePoint);
    for i = 1:num_frePoint
        pB_d = HB_Desired(:, i); % M_B * 1
        % closed solution
        % pD_d = zeros(size(HD_ctrl, 1), 1); % M_D*1
        % pd = [pB_d; pD_d];
        HB = squeeze(HB_ctrl(:, :, i));HD = squeeze(HD_ctrl(:, :, i));
        RB =  HB' * HB; RD = HD'* HD;
        % compute Filter coefficient 
        HBHpBd = HB'*pB_d;
        belta = 1e-6*max(eig(RB + RD));
        % belta = 1e-3;  % max(eig(R))*(1e-5);
        w(:, i) = (RB + RD + belta*eye(L))\HBHpBd;
    end
end
