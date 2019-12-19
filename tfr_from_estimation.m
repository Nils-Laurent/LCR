function [TFR_denoised] = tfr_from_estimation(sigma_s, STFT, phipE, phippE, L, Nfft)

% Assuming we have the exact ridge
% Cs = round(phipE/(L/Nfft)) + 1;
TFR_denoised = zeros(Nfft, L);
for n = 1:L
    Zn = -pi*sigma_s^2*(1 + 1i*phippE(n)*sigma_s^2)...
        /(1 + phippE(n)^2*sigma_s^4);
    Wk = min(Nfft, round(phipE(n)*Nfft/L)+1);
    %[~, Wk] = max(abs(STFT(:, n)));
    for k = 0:Nfft-1
        TFR_denoised(k+1, n) = STFT(Wk, n)...
            *exp(Zn*(k*L/Nfft - round(phipE(n)))^2);
    end
end

end

