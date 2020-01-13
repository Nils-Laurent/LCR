close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

f2_amp = 1024;
phi1 = 1000*t+f2_amp*(t.^2)/2;
phi2 = 1065*t+f2_amp*(t.^2)/2;

NR_phip = [phi1 phi2];
NR_phipp = [f2_amp*ones(L, 1) f2_amp*ones(L, 1)];

s_clean = exp(2*1i*pi*phi1) + exp(2*1i*pi*phi2);
NRidges = 2;
clwin = 1;

%% STFT parameters
Nfft = 512;
sigma = 1/sqrt(f2_amp);

%% Test
SNRs = -10:5:30;
NRep = 10;

% stft test
[g, Lg] = create_gaussian_window(L, Nfft, sigma);
[TFR_clean] = tfrstft(s_clean, Nfft, 1, g, Lg);
figure;
imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_clean));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square
xlabel('time');
ylabel('frequency');
pause

[SNR_LC, SNR_SSR_HT_init, SNR_SSR_HT_var, SNR_HT] = test_denoise_LC(...
    s_clean, NRidges, clwin, NR_phip, NR_phipp, sigma, Nfft, SNRs, NRep);

figure;
hold on;
plot(SNRs, SNR_LC, 'k-o', 'DisplayName', 'LC');
plot(SNRs, SNR_SSR_HT_init, ':', 'DisplayName', 'SSR-HT-init');
plot(SNRs, SNR_SSR_HT_var, '--', 'DisplayName', 'SSR-HT-var');
plot(SNRs, SNR_HT, 'DisplayName', 'HT');
hold off;
xlabel('input SNR');
ylabel('output SNR');
legend;
savefig('fig_2lin_autosave');


