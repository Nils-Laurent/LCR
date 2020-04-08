close all;

%% signal definition
L = 4096;
t = (0:L-1)'/L;

phi_lin = 250*t+568*(t.^2)/2;
s_lin = exp(2*1i*pi*phi_lin);
sigma_lin = 1/sqrt(568);

phi_cos = 1200*t + 60*cos(2*pi*t);
s_cos = exp(2*1i*pi*phi_cos);
sigma_cos = 0.0267;

% B = log(4096 - 510); % IF < L when t = 1
% phi_exp = 500*t+exp(B*t)/B;
% s_exp = exp(2*1i*pi*phi_exp);
% sigma_exp = 0.028;

phi_cos2 = 2000*t + 238*cos(2*pi*t);
s_cos2 = exp(2*1i*pi*phi_cos2);
sigma_cos2 = 0.0134;

Nfft = 512;
SNR_IN = -5:5:30;
NRep = 10;
clwin = 10;

%% figure global settings
set(groot, 'defaultLegendInterpreter', 'latex');

%% LCR test : linear chrip
[SNR_LCR, SNR_HT, SNR_SSR_HT,...
 SNR_LCR_c2, SNR_HT_c2, SNR_SSR_HT_c2, SNR_SST2] = test_LCR(...
    s_lin, clwin, sigma_lin, Nfft, SNR_IN, NRep);

plot_SNR_ALL(SNR_IN, SNR_LCR, SNR_HT, SNR_SSR_HT,...
    SNR_LCR_c2, SNR_HT_c2, SNR_SSR_HT_c2, SNR_SST2);
savefig('autosave_fig2_lin');
saveas(gcf,'autosave_fig2_lin','epsc');
close all;

%% LCR test : cos
[SNR_LCR, SNR_HT, SNR_SSR_HT,...
 SNR_LCR_c2, SNR_HT_c2, SNR_SSR_HT_c2, SNR_SST2] = test_LCR(...
    s_cos, clwin, sigma_cos, Nfft, SNR_IN, NRep);

plot_SNR_ALL(SNR_IN, SNR_LCR, SNR_HT, SNR_SSR_HT,...
    SNR_LCR_c2, SNR_HT_c2, SNR_SSR_HT_c2, SNR_SST2);
savefig('autosave_fig2_cos');
saveas(gcf,'autosave_fig2_cos','epsc');
close all;

%% LCR test : cos2
[SNR_LCR, SNR_HT, SNR_SSR_HT,...
 SNR_LCR_c2, SNR_HT_c2, SNR_SSR_HT_c2, SNR_SST2] = test_LCR(...
    s_cos2, clwin, sigma_cos2, Nfft, SNR_IN, NRep);

plot_SNR_ALL(SNR_IN, SNR_LCR, SNR_HT, SNR_SSR_HT,...
    SNR_LCR_c2, SNR_HT_c2, SNR_SSR_HT_c2, SNR_SST2);
savefig('autosave_fig2_cos2');
saveas(gcf,'autosave_fig2_cos2','epsc');
close all;

%% plot function
function plot_SNR_ALL(SNR_IN, SNR_LCR, SNR_HT, SNR_SSR_HT,...
    SNR_LCR_c2, SNR_HT_c2, SNR_SSR_HT_c2, SNR_SST2)
    figure;
    hold on;
    plot(SNR_IN, SNR_LCR, 'k-^',...
        'DisplayName', 'LCR $M_1$',...
        'MarkerSize', 10, 'LineWidth', 2);
    plot(SNR_IN, SNR_LCR_c2, 'k--^',...
         'DisplayName', 'LCR $M_2$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_SSR_HT, '-o', 'Color', [0 0.4470 0.7410],...
        'DisplayName', 'SSR-HT $M_1$',...
        'MarkerSize', 10, 'LineWidth', 2);
    plot(SNR_IN, SNR_SSR_HT_c2, '--o', 'Color', [0 0.4470 0.7410],...
        'DisplayName', 'SSR-HT $M_2$',...
        'MarkerSize', 10, 'LineWidth', 2);
    
    plot(SNR_IN, SNR_HT, '-s', 'Color', [0.6350 0.0780 0.1840],...
        'DisplayName', 'HT $M_1$',...
        'MarkerSize', 10, 'LineWidth', 2);
    plot(SNR_IN, SNR_HT_c2, '--s', 'Color', [0.6350 0.0780 0.1840],...
        'DisplayName', 'HT $M_2$',...
        'MarkerSize', 10, 'LineWidth', 2);

    plot(SNR_IN, SNR_SST2, '-x', 'Color', [0.4660 0.6740 0.1880],...
        'DisplayName', 'SST2',...
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
    set(gcf, 'Position',  [0, 0, 1000, 1000])
end