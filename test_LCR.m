function [SNR_LCR, SNR_HT, SNR_SSR_HT,...
    SNR_LCR_c2, SNR_HT_c2, SNR_SSR_HT_c2, SNR_SST2] = test_LCR(...
    modes, clwin, sigma_s, Nfft, SNR_IN, NRep)

Nr = size(modes, 2);
s_in = sum(modes, 2);
L = size(s_in, 1);

SNR_LCR = zeros(length(SNR_IN), Nr);
SNR_HT = zeros(length(SNR_IN), Nr);
SNR_SSR_HT = zeros(length(SNR_IN), Nr);

SNR_LCR_c2 = zeros(length(SNR_IN), Nr);
SNR_HT_c2 = zeros(length(SNR_IN), Nr);
SNR_SSR_HT_c2 = zeros(length(SNR_IN), Nr);

SNR_SST2 = zeros(length(SNR_IN), Nr);

for k=1:length(SNR_IN)
    for l=1:NRep
        fprintf("snr %d/%d, rep %d/%d\n", k, length(SNR_IN), l, NRep);
        
        noise = randn(L, 1)+1i*randn(L, 1);
        s_noise = sigmerge(s_in, noise, SNR_IN(k));
        
        fprintf("LCR-1, ");
        %% Model based denoising cas 1
        [m_LC, ~, Lg, ~] = denoise_LCR(...
            s_noise, Nr, clwin, sigma_s, Nfft, 1);
        
        fprintf("SSR-HT-1, ");
        %% SSR-HT denoising cas 1
        [m_HT, m_SSR_HT, ~] = denoise_SSR_HT_variant(...
            s_noise, Nr, clwin, Nfft, 'gauss', sigma_s, 1);
        
        fprintf("LCR-2, ");
        %% Model based denoising cas 2
        [m_LC_c2, ~, ~, ~] = denoise_LCR(...
            s_noise, Nr, clwin, sigma_s, Nfft, 2);
        
        fprintf("SSR-HT-2, ");
        %% SSR-HT denoising cas 2
        [m_HT_c2, m_SSR_HT_c2, ~] = denoise_SSR_HT_variant(...
            s_noise, Nr, clwin, Nfft, 'gauss', sigma_s, 2);
        
        fprintf("SST2\n");
        %% SST2
        [m_SST2, ~] = denoise_SST2(s_noise, Nr, clwin, sigma_s, Nfft);

        
        %% SNR signal
        X_win = 2*Lg:(L-2*Lg);
        for r = 1:Nr
            ref = modes(X_win, r);
            x_LCR = snr(ref, m_LC(X_win, r) - ref);
            x_HT = snr(ref, m_HT(X_win, r) - ref);
            x_SST_HT_var = snr(ref, m_SSR_HT(X_win, r) - ref);

            x_LCR_c2 = snr(ref, m_LC_c2(X_win, r) - ref);
            x_HT_c2 = snr(ref, m_HT_c2(X_win, r) - ref);
            x_SST_HT_var_c2 = snr(ref, m_SSR_HT_c2(X_win, r) - ref);

            x_SST2 = snr(ref, m_SST2(X_win, r) - ref);
        
            SNR_LCR(k, r) = SNR_LCR(k, r) + x_LCR;
            SNR_HT(k, r) = SNR_HT(k, r) + x_HT;
            SNR_SSR_HT(k, r) = SNR_SSR_HT(k, r) + x_SST_HT_var;

            SNR_LCR_c2(k, r) = SNR_LCR_c2(k, r) + x_LCR_c2;
            SNR_HT_c2(k, r) = SNR_HT_c2(k, r) + x_HT_c2;
            SNR_SSR_HT_c2(k, r) = SNR_SSR_HT_c2(k, r) + x_SST_HT_var_c2;

            SNR_SST2(k, r) = SNR_SST2(k, r) + x_SST2;
        end
    end
end

SNR_LCR = SNR_LCR/NRep;
SNR_HT = SNR_HT/NRep;
SNR_SSR_HT = SNR_SSR_HT/NRep;

SNR_LCR_c2 = SNR_LCR_c2/NRep;
SNR_HT_c2 = SNR_HT_c2/NRep;
SNR_SSR_HT_c2 = SNR_SSR_HT_c2/NRep;

SNR_SST2 = SNR_SST2/NRep;

end

