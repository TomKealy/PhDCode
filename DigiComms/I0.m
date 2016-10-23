function f=I0(beta)
%Bessel function of order 0
exp_bcos=inline('exp(beta*cos(theta))','theta','beta');
pi2=2*pi;
for n=1:length(beta)
   f(n)=quad('????????',0,pi2,0.0001,[],beta(n))/pi2;
end 
