function p=Hermitp(N)
%Hn+1(x)=2xHn(x)-Hn'(x) from 'Advanced Engineering Math' by Kreyszig     
if N<=0,  p=1;
 else 
   p=[2 0];
   for n=2:N, p= 2*[p 0]-[0 0 polyder(p)]; end
end
