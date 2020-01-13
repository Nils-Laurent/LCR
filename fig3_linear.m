close all

%% signal definition
L = 4096;
t = (0:L-1)'/L;

CRs = [568, 1200, 3500];
%CRs = 1200;
f1_amp = 250;
clwin = 10;

%% STFT parameters
Nfft = 512;

%% Test
SNRs = -10:5:30;
NRep = 10;

TFR_all = zeros(Nfft, L);
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

    % stft test
    [g, Lg] = create_gaussian_window(L, Nfft, sigma);
    [TFR_sig] = tfrstft(s_clean, Nfft, 1, g, Lg);
    TFR_all = TFR_all + TFR_sig;
end

figure;
imagesc((0:L-1)/L, (L/Nfft)*(1:Nfft), abs(TFR_all));
set(gca,'ydir','normal');
colormap(flipud(gray));
axis square
xlabel('time');
ylabel('frequency');
pause

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
    file_name = sprintf("lin_%d_autosave", CR);
    savefig(file_name);
%     pause
end
