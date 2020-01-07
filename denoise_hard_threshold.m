function [s_denoised] = denoise_hard_threshold(s_noise, sigma, Nfft)

cas = 1;
L = length(s_noise);
[g, Lg] = create_gaussian_window(L, Nfft, sigma);
[TFR_noise] = tfrstft(s_noise, Nfft, cas, g, Lg);
TFR_noise = TFR_noise/L;

Y2 = real(TFR_noise);
gamma = median(abs(Y2(:)))/0.6745;
TFR_th = TFR_noise;
TFR_th(abs(TFR_th) < 3*gamma) = 0;

[s_denoised] = L*itfrstft(TFR_th, cas, g);

end

