function [out_SNR_LC, out_SNR_SSR_HT_init, out_SNR_SSR_HT_var, out_SNR_HT] = test_denoise_LC(...
    signal, NRidges, clwin, NR_phip, NR_phipp, sigma, Nfft, SNRs, NRep)
L = size(signal, 1);

out_SNR_LC = zeros(length(SNRs), 1);
out_SNR_SSR_HT_init = zeros(length(SNRs), 1);
out_SNR_SSR_HT_var = zeros(length(SNRs), 1);
out_SNR_HT = zeros(length(SNRs), 1);

for k=1:length(SNRs)
    for l=1:NRep
        fprintf("snr %d/%d, rep %d/%d\n", k, length(SNRs), l, NRep);
        
        WGN = randn(L, 1)+1i*randn(L, 1);
        x = sigmerge(signal, WGN, SNRs(k));
        
        %% Model based denoising
        [s_LC, ~, Lg, ~] = denoise_LC(x, NRidges, clwin, sigma, Nfft, NR_phip, NR_phipp);
        
        %% SSR-HT denoising
        [s_HT, s_SSR_HT_var, s_SSR_HT_init] = denoise_SSR_HT_variant(x, NRidges, clwin, Nfft, 'gauss', sigma);
        
        %% SNR signal
        X_win = 2*Lg:(L-2*Lg);
        SNR_LC = snr(signal(X_win), s_LC(X_win) - signal(X_win));
        SNR_SSR_HT_init = snr(signal(X_win), s_SSR_HT_init(X_win) - signal(X_win));
        SNR_SSR_HT_var = snr(signal(X_win), s_SSR_HT_var(X_win) - signal(X_win));
        SNR_TH = snr(signal(X_win), s_HT(X_win) - signal(X_win));
        
        out_SNR_LC(k) = out_SNR_LC(k) + SNR_LC;
        out_SNR_SSR_HT_init(k) = out_SNR_SSR_HT_init(k) + SNR_SSR_HT_init;
        out_SNR_SSR_HT_var(k) = out_SNR_SSR_HT_var(k) + SNR_SSR_HT_var;
        out_SNR_HT(k) = out_SNR_HT(k) + SNR_TH;
    end
end

out_SNR_LC = out_SNR_LC/NRep;
out_SNR_SSR_HT_init = out_SNR_SSR_HT_init/NRep;
out_SNR_SSR_HT_var = out_SNR_SSR_HT_var/NRep;
out_SNR_HT = out_SNR_HT/NRep;

end

