function test_denoise_mode_est(s_noise, s_clean, phi, phip, phipp, amp)

L = length(s_noise);

Nfft = L;
cas = 1;

sigma_t = 0.0316*L/Nfft;
[g, Lg] = create_gaussian_window(Nfft, sigma_t);
sigma_s = sigma_t*Nfft/L;
% plot(1:2*Lg+1,g,1:2*Lg+1,exp(-pi*((-Lg:Lg)/N).^2/(sigma_s*Nfft/N)^2),'--');
% pause

% [TFR_clean] = tfrstft(s_clean, Nfft, cas, g, Lg);
% TFR_clean = TFR_clean/L;

[TFR_noise] = tfrstft(s_noise, Nfft, cas, g, Lg);
TFR_noise = TFR_noise/L;

% figure;
% imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_noise));
% set(gca,'ydir','normal');
% axis square
% colorbar;
% % return;


%% ridge extraction
[Cs] = exridge(TFR_noise,0,0,50);

% figure;
% hold on;
% imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_noise));
% set(gca,'ydir','normal');
% axis square
% colorbar;
% plot((0:L-1)/L, (L/Nfft)*Cs,'r','linewidth',1);
% hold off;
% return;

% Assuming we have the exact ridge
% Cs = round(phip/(L/Nfft)) + 1;

%% mode estimation

% TFR_backg = zeros(size(TFR_clean));
[~, ampE, phiE, phipE, phippE, ~, ~] = retrieve_mode(s_noise, Nfft, g, Lg, sigma_s, Cs);

%% plot estimates
% range_X = (Lg:L-Lg);

% figure;
% hold on;
% plot(range_X-1, ampE(range_X), 'DisplayName', "AE");
% plot(range_X-1, amp(range_X), 'DisplayName', "A'");
% legend;
% hold off;
% pause

% figure;
% hold on;
% plot(range_X-1, phipE(range_X), 'DisplayName', "phi'E");
% plot(range_X-1, phip(range_X), 'DisplayName', "phi'");
% plot(range_X-1, (Cs(range_X)-1)*L/Nfft, 'DisplayName', "Cs");
% legend;
% hold off;

% figure;
% hold on;
% plot(range_X-1, phippE(range_X), 'DisplayName', "phi''E");
% plot(range_X-1, phipp(range_X), 'DisplayName', "phi''");
% legend;
% hold off;

%% use estimate and inverse STFT

X = 2*Lg:(L-2*Lg);

[TFR_denoised] = tfr_from_estimation(sigma_s, TFR_noise, phipE, phippE, L, Nfft);
[s_denoised] = itfrstft(TFR_denoised, cas, g);
s_denoised = L*s_denoised;
out_snr_estimation = snr(s_clean(X), s_denoised(X) - s_clean(X));

%% Threshold denoising
TFR_th = TFR_noise;
Y2 = real(TFR_noise);
gamma = median(abs(Y2(:)))/0.6745;
TFR_th(abs(TFR_th) < 3*gamma) = 0;

[s_denoised_th] = itfrstft(TFR_th, cas, g);
s_denoised_th = L*s_denoised_th;
out_snr_th = snr(s_clean(X), s_denoised_th(X) - s_clean(X));

%% slice STFT

% figure;
% subplot(2, 1, 1);
% hold on;
% plot((0:Nfft-1)*L/Nfft,real(TFR_denoised(:,L/2)), 'DisplayName', "denoised");
% plot((0:Nfft-1)*L/Nfft,real(TFR_clean(:,L/2)),'--', 'DisplayName', "clean");
% hold off;
% subplot(2, 1, 2);
% hold on;
% plot((0:Nfft-1)*L/Nfft,imag(TFR_denoised(:,L/2)), 'DisplayName', "denoised");
% plot((0:Nfft-1)*L/Nfft,imag(TFR_clean(:,L/2)),'--', 'DisplayName', "clean");
% legend;
% hold off;

%% display

fprintf('______________________________\n');
fprintf('none   = %d \n', out_snr_estimation);
fprintf('------\n');
fprintf('TH     = %d \n', out_snr_th);

end
