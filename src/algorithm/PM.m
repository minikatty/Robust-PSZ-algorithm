function[w] = PM(Hb1, Hd1, HbDesired1, para)
%     freq = para.freq;
    frePoint = para.frePoint;
    for i = frePoint
        pdB = HbDesired1(:, i);
        % closed solution
        pdD = zeros(size(Hd1, 1), 1);
        pd = [pdB; pdD];
        Hbd = [squeeze(Hb1(:, :, i)); squeeze(Hd1(:, :, i))];
        % compute Filter coefficient 
        R = (Hbd') * Hbd;
        belta = 0;%max(eig(R))*(1e-5);
        w(:, i) = inv(R+belta*eye(size(Hb1, 2))) * (Hbd') * pd;
    end
end
%%
% function[w] = PM(Hb1, Hd1, HbDesired1)
%     pdB = HbDesired1;
%     % closed solution
%     pdD = zeros(size(Hd1, 1), 1);
%     pd = [pdB; pdD];
%     Hbd = [Hb1(:, :); Hd1(:, :)];
%     % compute Filter coefficient with Inevitable solution method
%     R = (Hbd') * Hbd;
%     belta = 0;%max(eig(R))*(1e-5);
%     w(:, 1) = inv(R+belta*eye(size(Hb1, 2))) * (Hbd') * pd;
% 
% end
% 
