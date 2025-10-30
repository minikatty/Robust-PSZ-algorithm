function[w] = ACC(Hb_ctrl, Hd_ctrl, target_indices)
    % Mb = size(Hb_ctrl, 1);
    % Md = size(Hd_ctrl, 1);
    % Mb and Md 始终是与w无关的常数
    L  = size(Hb_ctrl, 2);
    w = zeros(L,length(target_indices)); % [Loudspeaker,f]
    for i = 1:length(target_indices)
        idx = target_indices(i);
        % compute Filter coefficient with Inevitable solution method
        Rb = (squeeze(Hb_ctrl(:, :, idx))') * squeeze(Hb_ctrl(:, :, idx));
        Rd = (squeeze(Hd_ctrl(:, :, idx))') * squeeze(Hd_ctrl(:, :, idx));
        belta = max(eig(Rd))*(1e-8);
%         belta = max(eig(Rd))*(1e-5);
        R = (Rd + belta*eye(size(Hd_ctrl, 2))) \ (Rb);
        [~, w(:, i)] = MaxEigenvector(R);
    end
    
end

