function [out_SNR_phipE1, out_SNR_phipE2] = test_estimates_LC(...
    signal, NRidges, clwin, NR_phip, NR_phipp, sigma, Nfft, SNRs, NRep)
L = size(signal, 1);

out_SNR_phipE2 = zeros(length(SNRs), NRidges);
out_SNR_phipE1 = zeros(length(SNRs), NRidges);

for k=1:length(SNRs)
    for l=1:NRep
        fprintf("snr %d/%d, rep %d/%d\n", k, length(SNRs), l, NRep);
        
        WGN = randn(L, 1)+1i*randn(L, 1);
        x = sigmerge(signal, WGN, SNRs(k));
        
        %% Model based denoising
        [~, ~, Lg, E2] = denoise_LC(x, NRidges, clwin, sigma, Nfft, NR_phip, NR_phipp);
        
        %% SNR signal
        X_win = 2*Lg:(L-2*Lg);
        
        %% SNR parameters
        for r=1:NRidges
            C_phip = NR_phip(X_win, r);
            E2r = squeeze(E2(r, :, X_win));
            E2r = transpose(E2r);
            out_SNR_phipE2(k, r) = out_SNR_phipE2(k, r) ...
                + snr(C_phip, E2r(:, 2) - C_phip);
            out_SNR_phipE1(k, r) = out_SNR_phipE1(k, r) ...
                + snr(C_phip,E2r(:, 1) - C_phip);
        end
    end
end

out_SNR_phipE2 = out_SNR_phipE2/NRep;
out_SNR_phipE1 = out_SNR_phipE1/NRep;

end

