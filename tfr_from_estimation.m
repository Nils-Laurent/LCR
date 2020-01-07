function [TFR_denoised] = tfr_from_estimation(sigma, gamma, STFT, phipE, phippE, L, Nfft)

TFR_denoised = zeros(Nfft, L);
for n = 1:L
    Zn = pi*sigma^2*(1 + 1i*phippE(n)*sigma^2)...
        /(1 + phippE(n)^2*sigma^4);
    
    KWd = min(Nfft, round(phipE(n)*Nfft/L)+1);
    %Wd = ((KWd-1)*L/Nfft);
    %Gn = STFT(KWd, n)*exp(Zn*(Wd - phipE(n))^2);
    
    %% set threshold depending on the window
    th_sigma = 1/sqrt(2*pi)*sqrt(1/sigma^2 + sigma^2*real(phippE(n))^2);
    th_sigma = round(th_sigma*Nfft/L);
    
    %% reduce threshold to exclude noise
    th = 0;
    for k = 1:th_sigma
        if (STFT(KWd + th_sigma, n) < 3*gamma) ||...
                (STFT(KWd - th_sigma, n) < 3*gamma)
            break;
        end
        th = k;
    end
    
    %% compute |stft(n, phi'(n))|
    lower = max(1, KWd - th);
    upper = min(Nfft, KWd + th);
    Gn = mean(STFT(lower:upper, n)...
        .*exp(Zn*((lower-1:upper-1)'*L/Nfft - phipE(n)).^2));
    
    %% assign coefficients
    for k = 1:Nfft
        TFR_denoised(k, n) = Gn*exp(-Zn*((k-1)*L/Nfft - phipE(n))^2);
    end
end

end

