function y=QAM4_slicer(x)
% To find one of 4-QAM symbols which is the closest to x
PAM2= [-1  1]/sqrt(2); % 2 PAM levels
y= zeros(size(x));
for n=1:length(x)
   [min1,i1]= min(abs(PAM2-real(x(n))));
   [min2,i2]= min(abs(PAM2-imag(x(n))));
   y(n)= PAM2(i1) + j*PAM2(i2);
end