close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

f2_amp = 1024;
phi1 = 1000*t+f2_amp*(t.^2)/2;
phi2 = 1065*t+f2_amp*(t.^2)/2;

NR_phip = [phi1 phi2];
NR_phipp = [f2_amp*ones(L, 1) f2_amp*ones(L, 1)];

s1 = exp(2*1i*pi*phi1);
s2 = exp(2*1i*pi*phi2);
s_clean = s1 + s2;
NRidges = 2;
clwin = 1;

%% STFT parameters
Nfft = 512;
sigma = 1/sqrt(f2_amp);

%% Test
SNRs = -10:5:30;
NRep = 10;

all_l2 = zeros(length(SNRs), 1);
for k=1:length(SNRs)
    for l=1:NRep
        fprintf("snr %d/%d, rep %d/%d\n", k, length(SNRs), l, NRep);
        
        WGN = randn(L, 1)+1i*randn(L, 1);
        x = sigmerge(s_clean, WGN, SNRs(k));
        
        %% SSR-HT denoising
        [~, ~, ~, modes] = denoise_SSR_HT_variant(x, NRidges, clwin, Nfft, 'gauss', sigma);
        
        %% SNR signal
        X_win = 2*Lg:(L-2*Lg);
        l2_modes = norm(modes(X_win, 1) - modes(X_win, 2), 2);
        
        all_l2(k) = all_l2(k) + l2_modes;
    end
end
all_l2 = all_l2/NRep;

figure;
plot(SNRs, all_l2, 'DisplayName', '$\| \tilde{f_1} - \tilde{f_2} \|_2$');
xlabel('input SNR');
ylabel('l2 norm');
legend;
savefig('fig_2lin_diff_autosave');


