function  [c,se]=ade(c,y,a,delta)
%adaptive equalizer
%input:  y=the output response of channel 
%           a=message data
%     delta=the step size 
%output: c=updated equalizer coefficients
e=c*y'-a; se=e^2;
c=c-delta*e*y;
