function [modes, TFR_denoised, Lg, E2] = denoise_LCR(s_noise, Nr, clwin, sigma_s, Nfft, cas)

L = length(s_noise);
% cas = 1;

[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);


[TFR_noise] = tfrstft(s_noise, Nfft, cas, g, Lg);
TFR_noise = TFR_noise/L;

%% ridge extraction
[Cs] = exridge_mult(TFR_noise, Nr, 0, 0, clwin);

TFR_denoised = zeros(size(TFR_noise));
E2 = zeros(Nr, 2, L);

modes = zeros(L, Nr);

for r = 1:Nr
    %% mode estimation
    [phipE1, phipE2, phippE] = compute_FM(s_noise, Nfft, g, Lg, sigma_s, Cs(r, :));
    X = [phipE1, phipE2];
    E2(r, :, :) = transpose(X);

    %% use estimate and inverse STFT
    [TFR_denoised_r] = tfr_from_estimation(sigma_s, TFR_noise, phipE2, phippE, L, Nfft);
    TFR_denoised = TFR_denoised + TFR_denoised_r;
    modes(:, r) = L*itfrstft(TFR_denoised_r, cas, g);
end

end
