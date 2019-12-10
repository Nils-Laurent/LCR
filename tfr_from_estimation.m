function [TFR_denoised] = tfr_from_estimation(sigma_s, STFT, phipE, phippE, L, Nfft)

TFR_denoised = size(Nfft, L);
for n = 1:L
    R_n = -pi*sigma_s^2*(1 + 1i*phippE(n)*sigma_s^2)...
        /(1 + phippE(n)^2*sigma_s^4);
    for k = 0:Nfft-1
        TFR_denoised(k+1, n) = STFT(round(phipE(n)*Nfft/L) + 1, n)...
            *exp(R_n*(k*L/Nfft - phipE(n))^2);
    end
end

end

