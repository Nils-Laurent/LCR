close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

CR = 1024;
phi1 = 1000*t+CR*(t.^2)/2;
phi2 = 1065*t+CR*(t.^2)/2;
s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
modes = zeros(L, 2);
modes(:, 1) = s1;
modes(:, 2) = s2;
s_clean = s1 + s2;
NRidges = 2;
clwin = 1;

%% STFT parameters
Nfft = 512;
sigma_s = 1/sqrt(CR);

%% TFR
[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);
[TFR_clean] = tfrstft(s_clean, Nfft, 1, g, Lg);

figure;
imagesc((0:L-1)/L, (L/Nfft)*(0:Nfft-1), abs(TFR_clean));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square
xlabel('time', 'interpreter', 'latex');
ylabel('frequency', 'interpreter', 'latex');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
pause

%% MR test with mode mixing
SNR_IN = -5:5:30;
NRep = 10;
[SNR_LCR, SNR_HT, SNR_SSR_HT,...
 SNR_LCR_c2, SNR_HT_c2, SNR_SSR_HT_c2, SNR_SST2] = test_LCR(...
    modes, clwin, sigma_s, Nfft, SNR_IN, NRep);

set(groot, 'defaultLegendInterpreter','latex');
figure;
hold on;

plot(SNR_IN, SNR_LCR(:, 1), 'k-^',...
    'DisplayName', 'LCR $M_1 f_1$',...
    'MarkerSize', 10, 'LineWidth', 2);
plot(SNR_IN, SNR_SSR_HT(:, 1), '-o', 'Color', [0 0.4470 0.7410],...
    'DisplayName', 'SSR-HT $M_1 f_1$',...
    'MarkerSize', 10, 'LineWidth', 2);
plot(SNR_IN, SNR_HT(:, 1), '-s', 'Color', [0.6350 0.0780 0.1840],...
    'DisplayName', 'HT $M_1 f_1$',...
    'MarkerSize', 10, 'LineWidth', 2);

plot(SNR_IN, SNR_SST2(:, 1), '-x', 'Color', [0.4660 0.6740 0.1880],...
    'DisplayName', 'SST2 $f_1$',...
    'MarkerSize', 10, 'LineWidth', 2);

plot(SNR_IN, SNR_LCR(:, 2), 'k-.^',...
    'DisplayName', 'LCR $M_1 f_2$',...
    'MarkerSize', 10, 'LineWidth', 2);
plot(SNR_IN, SNR_SSR_HT(:, 2), '-.o', 'Color', [0 0.4470 0.7410],...
    'DisplayName', 'SSR-HT $M_1 f_2$',...
    'MarkerSize', 10, 'LineWidth', 2);
plot(SNR_IN, SNR_HT(:, 2), '-.s', 'Color', [0.6350 0.0780 0.1840],...
    'DisplayName', 'HT $M_1 f_2$',...
    'MarkerSize', 10, 'LineWidth', 2);

plot(SNR_IN, SNR_SST2(:, 2), '-.x', 'Color', [0.4660 0.6740 0.1880],...
    'DisplayName', 'SST2 $f_2$',...
    'MarkerSize', 10, 'LineWidth', 2);
hold off;

xlabel('input SNR', 'interpreter', 'latex');
ylabel('output SNR', 'interpreter', 'latex');
xlim([SNR_IN(1), SNR_IN(length(SNR_IN))]);

lgd = legend('Location', 'northwest');
lgd.FontSize = 24;
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 26);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 26);
pbaspect([1 1 1]);
set(gcf, 'Position',  [0, 0, 1000, 1000]);
savefig('autosave_fig3_MR_LC');
saveas(gcf,'autosave_fig3_MR_LC','epsc');


