function[w] = ACC(HB_ctrl, HD_ctrl)
    % Mb = size(Hb_ctrl, 1);
    % Md = size(Hd_ctrl, 1);
    % Mb and Md 始终是与w无关的常数 
    [~, L, num_frePoint]  = size(HB_ctrl);
    w = zeros(L,num_frePoint); % [Loudspeaker,num_frePoint]
    for i = 1:num_frePoint
        % compute Filter coefficient with Inevitable solution method
        Rb = (squeeze(HB_ctrl(:, :, i))') * squeeze(HB_ctrl(:, :, i));
        Rd = (squeeze(HD_ctrl(:, :, i))') * squeeze(HD_ctrl(:, :, i));
        % belta = 0;
        % belta = max(eig(Rd))*1e-5;
        belta = 1e-9; % 
        R = (Rd + belta*max(eig(Rd))*eye(size(HD_ctrl, 2))) \ (Rb); % 
        % R = Rd\Rb;
        [~, w(:, i)] = MaxEigenvector(R);
    end 
end

