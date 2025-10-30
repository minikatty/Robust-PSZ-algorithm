function[w] = ACC_PM(Hb1, Hd1, pdB1, kappa, para) 
    frePoint = para.frePoint;
%     freq = para.freq;
    for i = frePoint
        Hb = squeeze(Hb1(:, :, i));
        Hd = squeeze(Hd1(:, :, i));
        pdB = pdB1(:, i);
        % compute Filter coefficient 
        R = kappa*(Hd')*Hd + (1-kappa)*(Hb')*Hb;
%         belta = 0;%max(eig(R))*(1e-5);
        belta = max(eig(R))*(2e-6);
        w(:, i) = inv(R + belta*eye(size(Hb, 2))) * (1-kappa) * (Hb') * pdB;  
    end
end
