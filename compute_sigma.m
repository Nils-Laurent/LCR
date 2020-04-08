%% signal definition
L = 4096;
t = (0:L-1)'/L;

phi_cos2 = 2000*t + 238*cos(2*pi*t);
s_clean = exp(2*1i*pi*phi_cos2);

Nfft = 4096;

alpha = 3;
sigma_set = 0.0132:0.0001:0.0136;
SL = length(sigma_set);
RE_vec = zeros(1, SL);
iSL = 0;
for sigma = sigma_set
    iSL = iSL + 1;
    fprintf("%u/%u\n", iSL, SL);
    [g, Lg] = create_gaussian_window(L, Nfft, sigma);
    [TFR, omega, omega2, q] = FM_operators(s_clean, Nfft, g, Lg, sigma);
    Y = abs(TFR);
    TFR_MS = sum(Y(:));
    RE_vec(iSL) = 1/(1 - alpha)*log2(sum(Y(:).^alpha)/(TFR_MS^alpha)) - log2(L);
end

figure;
plot(sigma_set, RE_vec);