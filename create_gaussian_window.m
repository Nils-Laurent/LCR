function [h, Lh] = create_gaussian_window(Nfft, sigma_w)

prec = 10^(-3);
L =  Nfft*sigma_w;
Lh = floor(L*sqrt(-log(prec)/pi))+1;
h = amgauss(2*Lh+1,Lh+1,L);

end