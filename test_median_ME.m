function test_median_ME(signal, NRidges, clwin, NR_phip, NR_phipp, sigma, Nfft)
L = size(signal, 1);
%SNRs = -10:5:30;
SNRs = -10:10:10;
NRep = 10;

out_SNR_phipER = zeros(length(SNRs), NRidges);
out_SNR_phipE1 = zeros(length(SNRs), NRidges);
out_SNR_phipE2 = zeros(length(SNRs), NRidges);
out_SNR_phippE = zeros(length(SNRs), NRidges);
out_SNR_phipEM = zeros(length(SNRs), NRidges);
out_SNR_phippEM = zeros(length(SNRs), NRidges);

for k=1:length(SNRs)
    for l=1:NRep
        fprintf("snr %d/%d, rep %d/%d\n", k, length(SNRs), l, NRep);
        
        WGN = randn(L, 1)+1i*randn(L, 1);
        x = sigmerge(signal, WGN, SNRs(k));
        
        %% Model based denoising
        [~, ~, Lg, E6] = denoise_model_estim(x, NRidges, clwin, sigma, Nfft, NR_phip, NR_phipp);
        
        %% SNR signal
        X_win = 2*Lg:(L-2*Lg);
        
        %% SNR parameters
        for r=1:NRidges
            C_phip = NR_phip(X_win, r);
            C_phipp = NR_phipp(X_win, r);
            E6r = squeeze(E6(r, :, X_win));
            E6r = transpose(E6r);
            out_SNR_phipER(k, r) = out_SNR_phipER(k, r) ...
                + snr(C_phip,E6r(:, 1) - C_phip);
            out_SNR_phipE1(k, r) = out_SNR_phipE1(k, r) ...
                + snr(C_phip,E6r(:, 2) - C_phip);
            out_SNR_phipE2(k, r) = out_SNR_phipE2(k, r) ...
                + snr(C_phip, E6r(:, 3) - C_phip);
            out_SNR_phippE(k, r) = out_SNR_phippE(k, r) ...
                + snr(C_phipp, E6r(:, 4) - C_phipp);
            out_SNR_phipEM(k, r) = out_SNR_phipEM(k, r) ...
                + snr(C_phip, E6r(:, 5) - C_phip);
            out_SNR_phippEM(k, r) = out_SNR_phippEM(k, r) ...
                + snr(C_phipp, E6r(:, 6) - C_phipp);
        end
    end
end

fprintf("generating figures\n", k, length(SNRs), l, NRep);

out_SNR_phipER = out_SNR_phipER/NRep;
out_SNR_phipE1 = out_SNR_phipE1/NRep;
out_SNR_phipE2 = out_SNR_phipE2/NRep;
out_SNR_phippE = out_SNR_phippE/NRep;
out_SNR_phipEM = out_SNR_phipEM/NRep;
out_SNR_phippEM = out_SNR_phippEM/NRep;

set(groot, 'defaultLegendInterpreter','latex');

for r=1:NRidges
    figure;
    
    subplot(2, 1, 1);
    title(sprintf("$\\phi_%d'$ estimations", r), 'Interpreter', 'latex');
    hold on;
    plot(SNRs, out_SNR_phipEM(:, r), 'k-o', 'DisplayName', "SNR-$\omega^{(2)}_M(t, r(t))$");
    plot(SNRs, out_SNR_phipE2(:, r), '--', 'DisplayName', "SNR-$\omega^{(2)}(t, r(t))$");
    plot(SNRs, out_SNR_phipE1(:, r), ':', 'DisplayName', "SNR-$\omega(t, r(t))$");
    plot(SNRs, out_SNR_phipER(:, r), 'DisplayName', "SNR-$r(t)$");
    hold off;
    legend;
    
    subplot(2, 1, 2);
    title(sprintf("$\\phi_%d''$ estimations", r), 'Interpreter', 'latex');
    hold on;
    plot(SNRs, out_SNR_phippEM(:, r), 'k-o', 'DisplayName', 'SNR-$Re(\widetilde{q}(t, r(t)))_M$');
    plot(SNRs, out_SNR_phippE(:, r), '--', 'DisplayName', 'SNR-$Re(\widetilde{q}(t, r(t)))$');
    hold off;
    legend;
end

end

