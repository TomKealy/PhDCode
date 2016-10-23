%dc01e01.m
%plots the CTFS spectra of rectangular/triangular waves
clear, clf
global P D
N=10; %Highest order of CTFS coefficient and representations
D=1; P=2; CTFS('rD_wave',P,N/2,221);
w0=2*pi/P; k1=linspace(-N/2,N/2); 
RD1=sinc(k1*w0*D/2/pi); %Spectrum
hold on, plot(k1,abs(RD1),':') %Envelope for the spectrum
axis([-N/2 N/2 -0.2 1.2])
P=2; CTFS('tri_wave',P,N/2,223);
w0=2*pi/P; 
Tri=sinc(k1*w0*D/2/pi); Tri1=Tri.*Tri; %Spectrum
hold on, plot(k1,Tri1,':') %Envelope for the spectrum
axis([-N/2 N/2 -0.2 1.2]) 
