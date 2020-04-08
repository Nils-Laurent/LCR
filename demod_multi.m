function [sp1_s,sp2_s,integ1,integ2] = demod_multi(s,SST,VSST,omega,omega2,t,nr,jump)   

%Computation of the ridges with no regularization, effect of zero padding 
 beta   = 0;
 lambda = 0;
 
 %computation of the different ridges
 [Cs2, Es] = exridge_mult(SST,nr,lambda,beta,jump);
 [Cs21,Es] = exridge_mult(VSST,nr,lambda,beta,jump);

 integ1 = zeros(size(Cs2));
 integ2 = zeros(size(Cs21));
 
 sp1_s  =  zeros(nr,length(s));
 sp2_s  =  zeros(nr,length(s));
 
 for k = 1:nr
 
  Y = [];
  YY = [];
  for kk = 1: length(t),
   Y  = [Y omega(Cs2(k,kk),kk)];
   YY  = [YY omega2(Cs21(k,kk),kk)];
  end
  
  integ1(k,:) = cumtrapz(t,Y);
  sp1_s(k,:)  = s.*exp(-2*1i*pi*(integ1(k,:)-104.*t));
  
  integ2(k,:) = cumtrapz(t,YY);
  sp2_s(k,:)  = s.*exp(-2*1i*pi*(integ2(k,:)-104.*t));
 
 end 
 