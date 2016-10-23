%dc02p01.m 
clear, clf 
wc=pi/2;  th0=pi/4; 
t=0:0.001:10; Nt=length(t); 
x=0:0.1:5; NB=50; B=5;
for m=1:2 
   A=2*(m-1); Ac=A*cos(th0); As=A*sin(th0); 
   for n=1:Nt 
      rc= Ac+randn; rs=As+randn; 
      r(n)= rc*cos(wc*t(n))-rs*sin(wc*t(n)); 
      z(n)=sqrt(rc^2+rs^2);   
   end 
   z_below_B=z(find(z<B)); Nz=length(z); 
   subplot(210+m)
   [ns,cs]=hist(z_below_B,NB); dx=cs(2)-cs(1); 
   f= ????????(x,A,1); 
   plot(x,f*dx,'r',cs,ns/Nz,':') 
end
