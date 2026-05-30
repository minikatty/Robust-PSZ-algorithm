Y = fft(squeeze(RIR_tmp(1,1,:)));
L = length(Y);
Y = abs(Y/L);
P1 = Y(1:floor(L/2)+1); % 统一截取单边频谱
f = fs*(0:floor(L/2))/L; % 统一构建频率向量
figure;
plot(f,20*log10(P1));