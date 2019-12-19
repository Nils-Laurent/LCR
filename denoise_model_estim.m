function [s_denoised_ME, TFR_denoised, Lg, E6] = denoise_model_estim(s_noise, NRidges, sigma_s, Nfft, phip, phipp)

L = length(s_noise);
cas = 1;

[g, Lg] = create_gaussian_window(L, Nfft, sigma_s);
%sigma_s = sigma_t*Nfft/L;


[TFR_noise] = tfrstft(s_noise, Nfft, cas, g, Lg);
TFR_noise = TFR_noise/L;

%% ridge extraction
[Cs] = exridge_mult(TFR_noise, NRidges, 0, 0, 10);

% figure;
% hold on;
% imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_noise));
% set(gca,'ydir','normal');
% axis square
% colorbar;
% plot((0:L-1)/L, (L/Nfft)*Cs,'r','linewidth',1);
% hold off;
% pause;

TFR_denoised = zeros(size(TFR_noise));
E6 = zeros(NRidges, 6, L);
for r = 1:NRidges
    %% mode estimation
    [phipE1, phipE2, phippE, phipEM, phippEM] = retrieve_mode(s_noise, Nfft, g, Lg, sigma_s, Cs(r, :));
    X = [L/Nfft*Cs(r, :)', phipE1, phipE2, phippE, phipEM, phippEM];
    E6(r, :, :) = transpose(X);

    %% use estimate and inverse STFT
    %[TFR_denoised_r] = tfr_from_estimation(sigma_s, TFR_noise, phipEM, phippEM, L, Nfft);
    [TFR_denoised_r] = tfr_from_estimation(sigma_s, TFR_noise, phip, phipp, L, Nfft);
    TFR_denoised = TFR_denoised + TFR_denoised_r;
end

% figure;
% imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_denoised));
% set(gca,'ydir','normal');
% axis square
% colorbar;
% pause;

% XLg = 2*Lg:L-2*Lg;
% TFR_D_mod = abs(TFR_noise(:, XLg)) - abs(TFR_denoised(:, XLg));
% TFR_D_real = real(TFR_noise(:, XLg)) - real(TFR_denoised(:, XLg));
% TFR_D_imag = imag(TFR_noise(:, XLg)) - imag(TFR_denoised(:, XLg));
% 
% TFR_D_mod_all = abs(TFR_noise) - abs(TFR_denoised);
% 
% figure;
% imagesc(abs(TFR_D_real) - abs(TFR_D_imag));
% title("|real error| - |imag error|");
% set(gca,'ydir','normal');
% colorbar;
% figure;
% imagesc(TFR_D_mod);
% title("modulus error");
% set(gca,'ydir','normal');
% colorbar;
% 
% XERR = zeros(L, 1);
% for index = 1:L
%     kRidge = round(phip(index)*Nfft/L)+1;
%     XERR(index) = TFR_D_mod_all(kRidge, index);
% end
% 
% figure;
% F_XERR = fft(XERR, Nfft);
% F_XERR = F_XERR(Nfft/2+1:Nfft);
% plot((0:Nfft/2-1)*(L/Nfft), abs(F_XERR));
% title("error frequency");
% [fval, freq] = max(F_XERR);
% freq = (freq-1)*L/Nfft;
% fprintf("max|fft(Xerr)| at %d Hz\n", freq);
% factor(freq)
% 
% T0 = 1770;
% figure;
% subplot(2, 1, 1);
% plot(1:Nfft,real(TFR_noise(:,T0)),'k', 1:Nfft,real(TFR_denoised(:,T0)),'b--');
% title(sprintf("real, T0 = %d", T0));
% subplot(2, 1, 2);
% plot(1:Nfft,imag(TFR_noise(:,T0)),'k', 1:Nfft,imag(TFR_denoised(:,T0)),'r--');
% title(sprintf("imag, T0 = %d", T0));
% figure;
% plot(1:Nfft,real(TFR_denoised(:,T0)) - real(TFR_noise(:,T0)), 1:Nfft,imag(TFR_denoised(:,T0)) - imag(TFR_noise(:,T0)), '--');
% title(sprintf("errors, T0 = %d", T0));
% % pause

[s_denoised_ME] = L*itfrstft(TFR_denoised, cas, g);

end
