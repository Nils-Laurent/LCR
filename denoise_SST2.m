function [modes_SST2, STFT] = denoise_SST2(s_noise, Nr, clwin, sigma_s, Nfft)

L = length(s_noise);
cas = 1;

[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);

[TFR_noise] = tfrstft(s_noise, Nfft, cas, g, Lg);
TFR_noise = TFR_noise/L;
gamma = median(abs(real(TFR_noise(:))))/0.6745;
aa = 1/sigma_s^2/L;
t = (0:L-1)/L;
% s_noise = s_noise(:);

modes_SST2 = zeros(L, Nr);
[STFT,SST,SST2,omega,omega2] = sst2_new(s_noise,aa,Nfft,gamma);

[Cs, Er] = exridge_mult(SST2, Nr, 0, 0, clwin);

for p=1:Nr
    modes_SST2(:, p) = (L/Nfft)*recmodes(SST2, Cs(p,:), 5);
end

end

