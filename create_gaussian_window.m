function [h, Lh] = create_gaussian_window(L, Nfft, sigma_w)

prec = 10^(-3);
Lw =  L*sigma_w;
Lh = floor(Lw*sqrt(-log(prec)/pi))+1;
h = amgauss(2*Lh+1,Lh+1,Lw);

if 2*Lh + 1 > Nfft
    disp("create_gaussian_window [Warning] : 2*Lh+1 > Nfft");
end

end