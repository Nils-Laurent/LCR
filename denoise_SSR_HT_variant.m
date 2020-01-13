function [s_noise_hard, s_noise_soft_var, s_noise_soft_init, modes_var, tfr_noise_soft_init, h,Lh] = denoise_SSR_HT_variant(sn, nr, clwin, Nfft, window, sigma_opt, s_clean)

%INPUT
%sig        : type of studied signal
%window     : window used 
%SNR        : SNR corresponding to the added noise

%OUTPUT
%tfr_noise_hard     : noisy STFT, hard thresholded
%tfr_noise_soft     : noisy STFT, soft thresholded

cas = 2;

%we build the filter h
if strcmp(window,'hamming')
    hlength=floor(Lg);%optimal window determined by Renyi entropy
    hlength=hlength+1-rem(hlength,2);%the length of the filter has to be odd
    h = tftb_window(hlength,window);
    [hrow,hcol]=size(h); 
    Lh=(hrow-1)/2;
else
    %the window is the Gaussian window    
%     prec = 10^(-3);
%     L =  sigma_opt*Nfft;
%     Lh = floor(L*sqrt(-log(prec)/pi))+1;
%     h = amgauss(2*Lh+1,Lh+1,L); 
    [h, Lh] = create_gaussian_window(length(sn), Nfft, sigma_opt);
end

[tfr_noise] = tfrstft(sn,Nfft,cas,h,Lh);
% [tfr] = tfrstft(s_clean,Nfft,cas,h,Lh);

Y2 = real(tfr_noise);
gamma_estime = median(abs(Y2(:)))/0.6745;
% imagesc(abs(tfr_noise));
% figure
% imagesc(abs(tfr_noise) > 3*gamma_estime);
% set(gca,'ydir','normal');
% pause
%clwin = 10;
[Cs] = exridge_mult(tfr_noise,nr,0,0,clwin);

Cs = Cs';

B = size(tfr_noise); 
Abstfr      = abs(tfr_noise);

s_noise_hard = zeros(size(sn));
s_noise_soft_init = zeros(size(sn));
s_noise_soft_var = zeros(size(sn));

modes_var = zeros(length(sn), nr);

%construction of the TF mask
for j=1:nr,
    tfr_noise_hard = zeros(B); 
    tfr_noise_soft_init = zeros(B);
    tfr_noise_soft_var = zeros(B);
    
    for r = 1:B(2)  
        val = 3*gamma_estime; %threshold for the transform depending on the noise level
        if (Abstfr(Cs(r,j),r) > val)
            k1 = 0;
            k2 = 0;
            eta1 = - 1;
            while (eta1 < 0)&&(Abstfr(Cs(r,j)-min(k1,Cs(r,j)-1),r) > val)
                if (k1 ~= Cs(r,j)-1)
                    k1 = k1+1;
                else
                    eta1 = k1;
                end
            end
            if (eta1 < 0)
                eta1 = k1-1;
            end
            eta2 = -1;
            while (eta2 < 0) && (Abstfr(Cs(r,j)+min(k2,B(1)-Cs(r,j)),r) > val)
                if (k2 ~= B(1)-Cs(r,j))
                    k2 = k2+1;   
                else
                    eta2 = k2;
                end
            end
            if (eta2 < 0)
                eta2 = k2;
            end
            eta = max(eta1,eta2); %we take the larger value
            indmin = max(1,Cs(r,j)-eta);
            indmax = min(B(1),Cs(r,j)+eta);
            
            tfr_noise_hard(indmin:indmax,r) = tfr_noise(indmin:indmax,r);  % hard thresholding 
            
            %Shift step     
            if (abs(tfr_noise(Cs(r,j),r))-abs(tfr_noise(Cs(r,j)-1,r)) >= 2) &&...
                    (abs(tfr_noise(Cs(r,j),r))-abs(tfr_noise(Cs(r,j)+1,r)) >= 2)...
                    && (Cs(r,j)-eta >= 1) && (Cs(r,j)+eta <= B(1))
                [Y,I] = max(abs(imag(tfr_noise(indmin:indmax,r))));
                shift = I - eta -1;
                tfr_noise_shift_imag = circshift(imag(tfr_noise(:,r)),shift);
                [Y,I] = max(abs(real(tfr_noise(indmin:indmax,r))));
                shift = I - eta -1;
                tfr_noise_shift_real = circshift(real(tfr_noise(:,r)),shift);
                tfr_noise_shift = tfr_noise_shift_real + 1i*tfr_noise_shift_imag;
                
                %Symmetry step
                tfr_noise_soft_init(indmin:indmax,r) = 1/2*(tfr_noise_shift(indmin:indmax)+tfr_noise_shift(indmax:-1:indmin));
                
                %Symmetry step for variant
                tfr_noise_soft_var(indmin:indmax,r) = 1/2*(tfr_noise(indmin:indmax,r)+tfr_noise(indmax:-1:indmin,r));
            else
                tfr_noise_shift = tfr_noise(:,r);
                
                %Symmetry step
                if abs(tfr_noise_shift(Cs(r,j)+1)) >= abs(tfr_noise_shift(Cs(r,j)-1))
                    for k = 0:eta
                        if ((Cs(r,j)-k) >= 1)&&((Cs(r,j)+k+1) <= B(1))
                            tfr_noise_soft_init(Cs(r,j)-k,r)= 1/2*(tfr_noise_shift(Cs(r,j)-k)+tfr_noise_shift(Cs(r,j)+k+1));
                            tfr_noise_soft_init(Cs(r,j)+k+1,r)= 1/2*(tfr_noise_shift(Cs(r,j)-k)+tfr_noise_shift(Cs(r,j)+k+1));
                        end
                    end
                else
                    for k = 0:eta
                        if ((Cs(r,j)-k-1) >= 1)&&((Cs(r,j)+k) <= B(1))   
                            tfr_noise_soft_init(Cs(r,j)-k-1,r)= 1/2*(tfr_noise_shift(Cs(r,j)-k-1)+tfr_noise_shift(Cs(r,j)+k));
                            tfr_noise_soft_init(Cs(r,j)+k,r)= 1/2*(tfr_noise_shift(Cs(r,j)-k-1)+tfr_noise_shift(Cs(r,j)+k));
                        end
                    end
                end
                tfr_noise_soft_var(:, r) = tfr_noise_soft_init(:, r);
            end
              
            %smoothing step
            X = [1:max(1,Cs(r,j)-eta-4) max(1,Cs(r,j)-eta):min(B(1),Cs(r,j)+eta) min(B(1),Cs(r,j)+eta+4):B(1)];
            X  = unique(X);
            Y  = tfr_noise_soft_init(X,r);
            YY = pchip(X,real(Y),1:B(1));
            ZZ = pchip(X,imag(Y),1:B(1));
            tfr_noise_soft_init(:,r) = YY' + 1i*ZZ';
            
            %smoothing variant
            Y  = tfr_noise_soft_var(X,r);
            YY = pchip(X,real(Y),1:B(1));
            ZZ = pchip(X,imag(Y),1:B(1));
            tfr_noise_soft_var(:,r) = YY' + 1i*ZZ';
        end
    end
%             figure
%             imagesc(abs(tfr_noise_soft));
%             figure
%             imagesc(abs(tfr));
%             figure
%             plot(1:B(1),real(tfr(:,500)),1:B(1),real(tfr_noise_soft(:,500)),'--', 1:B(1),real(tfr_noise_hard(:, 500)), '-o');
%             figure
%             plot(1:B(1),imag(tfr(:,500)),1:B(1),imag(tfr_noise_soft(:,500)),'--', 1:B(1),imag(tfr_noise_hard(:, 500)), '-o');
%             pause

    s_noise_hard = s_noise_hard + itfrstft(tfr_noise_hard,cas,h);
    s_noise_soft_init = s_noise_soft_init + itfrstft(tfr_noise_soft_init,cas,h);
    
    cmode = itfrstft(tfr_noise_soft_var,cas,h);
    s_noise_soft_var = s_noise_soft_var + cmode;
    modes_var(:, j) = cmode;
end

% figure;
% hold on;
% imagesc((0:N-1)/N, (N/Nfft)*(0:Nfft-1), abs(tfr_noise_soft));
% plot((0:N-1)/N, N/Nfft*Cs(:, 1));
% plot((0:N-1)/N, N/Nfft*Cs(:, 2));
% hold off;
% set(gca,'ydir','normal');
% title("TFR");
% axis square
% pause

% s_noise_hard = itfrstft(tfr_noise_hard,cas,h);
% s_noise_soft_init = itfrstft(tfr_noise_soft,cas,h);
% s_noise_soft_var = itfrstft(tfr_noise_soft_var,cas,h);
