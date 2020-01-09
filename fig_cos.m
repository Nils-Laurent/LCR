close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

A = 3*pi;
X = 1200;
Y = 60;
phi = X*t+Y*cos(A*t);

s_clean = exp(2*pi*1i*phi);
NRidges = 1;
NR_phip = X - Y*A*sin(A*t);
NR_phipp = - Y*A^2*cos(A*t);
clwin = 10;

%% STFT parameters
Nfft = 512;
%sigma = 0.15*Nfft/L;
%sigma = 0.00187;
sigma = 1/sqrt(A^2*Y);

%% Test
SNRs = -10:5:30;
NRep = 10;

%stft test
[g, Lg] = create_gaussian_window(L, Nfft, sigma);
[TFR] = tfrstft(s_clean, Nfft, 1, g, Lg);
figure;
imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square
xlabel('time');
ylabel('frequency');
pause

%% test denoise
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
savefig('fig_cos_autosave');

