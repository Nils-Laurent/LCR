close all;

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

[g, Lg] = create_gaussian_window(L, Nfft, sigma_lin);
[TFR_lin] = tfrstft(s_lin, Nfft, 1, g, Lg);

[g, Lg] = create_gaussian_window(L, Nfft, sigma_cos);
[TFR_cos] = tfrstft(s_cos, Nfft, 1, g, Lg);

[g, Lg] = create_gaussian_window(L, Nfft, sigma_cos2);
[TFR_cos2] = tfrstft(s_cos2, Nfft, 1, g, Lg);

figure;
imagesc((0:L-1)/L, (L/Nfft)*(0:Nfft-1), abs(TFR_lin));
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

figure;
imagesc((0:L-1)/L, (L/Nfft)*(0:Nfft-1), abs(TFR_cos));
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
set(gcf, 'Position',  [0, 0, 1000, 1000])

figure;
imagesc((0:L-1)/L, (L/Nfft)*(0:Nfft-1), abs(TFR_cos2));
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
set(gcf, 'Position',  [0, 0, 1000, 1000])