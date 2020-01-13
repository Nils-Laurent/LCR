close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

CRs = [568, 1200, 3500];
f1_amp = 250;
clwin = 10;

%% STFT parameters
Nfft = 512;

%% Test
SNRs = -10:5:30;
NRep = 10;

out_SNRs_o1 = zeros(length(CRs), length(SNRs));
out_SNRs_o2 = zeros(length(CRs), length(SNRs));

iCR = 0;
for CR=CRs
    iCR = iCR + 1;
    
    %% Signal definition
    phi = f1_amp*t+CR*(t.^2)/2;

    s_clean = exp(2*1i*pi*phi);
    NRidges = 1;
    NR_phip = f1_amp + CR*t;
    NR_phipp = CR*ones(L, 1);
    
    %% STFT related operations
    sigma = 1/sqrt(CR);

    %% denoising test
    ftitle = sprintf("Chirp rate = %d", CR);
    fprintf("%s\n", ftitle);
    [SNRs_phi_o1, SNRs_phi_o2] = test_estimates_LC(...
        s_clean, NRidges, clwin, NR_phip, NR_phipp, sigma, Nfft, SNRs, NRep);
    out_SNRs_o2(iCR, :) = SNRs_phi_o2;
    out_SNRs_o1(iCR, :) = SNRs_phi_o1;
end

set(groot, 'defaultLegendInterpreter','latex');

figure;
hold on;
plot(SNRs, out_SNRs_o2(3, :), 'r', 'DisplayName', "$f_1, \hat{\omega}_{\widetilde{f}}^{[2]}$");
plot(SNRs, out_SNRs_o1(3, :), 'r--', 'DisplayName', "$f_1, \hat{\omega}_{\widetilde{f}}$");
plot(SNRs, out_SNRs_o2(2, :), 'b-*', 'DisplayName', "$f_2, \hat{\omega}_{\widetilde{f}}^{[2]}$");
plot(SNRs, out_SNRs_o1(2, :), 'b--*', 'DisplayName', "$f_2, \hat{\omega}_{\widetilde{f}}$");
plot(SNRs, out_SNRs_o2(1, :), 'k-o', 'DisplayName', "$f_3, \hat{\omega}_{\widetilde{f}}^{[2]}$");
plot(SNRs, out_SNRs_o1(1, :), 'k--o', 'DisplayName', "$f_3, \hat{\omega}_{\widetilde{f}}$");
hold off;
xlabel('input SNR');
ylabel('output SNR');
legend;
savefig('fig_estimators_autosave');
