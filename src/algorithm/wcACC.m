function [wwcACC] = wcACC(HB, HD, Hbdesired, para)
    frePoint = para.frePoint;
    epsilonB = para.epsilonB;
    epsilonD = para.epsilonD;
    for i = frePoint
        HB1 = squeeze(HB(:, :, i));
        HD1 = squeeze(HD(:, :, i));
        Rb = HB1'*HB1;
        Rd = HD1'*HD1;
        gammaB = epsilonB(i)^2;
        gammaD = epsilonD(i)^2;
        A = inv(Rd + gammaD*eye(size(Rd))) * (Rb - gammaB*eye(size(Rb)));
        [~, wwcACC(:, i)] = MaxEigenvector(A);
    end
end

