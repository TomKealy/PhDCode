function d=detector_PSK(th,M)
%th=atan2(ycsk(2),ycsk(1));
if th<-pi/M, th=th+2*pi; end
[thmin,ind]=min(abs(th-2*pi/M*[0:M-1])); % Eq.(7.3.9) or Fig. 7.7 
d=ind-1;
