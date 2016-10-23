close all;
clear all;

reg1 =0;
reg2 =0;
reg3 = 0;

eta =sqrt(2)/2;
theta =2*pi*1/100;
Kp = [(4*eta*theta)/(1+2*eta*theta+theta^2)];
Ki = [(4*theta^2)/(1+2*eta*theta+theta^2)];
d_phi_1 = 1/20;
n_data = 100;

for nn =1:n_data
  phi1= reg1 + d_phi_1;
  phi1_reg(nn) = phi1;

  s1 =exp(j*2*pi*reg1);
  s2 =exp(j*2*pi*reg2);

  s1_reg(nn) =s1;
  s2_reg(nn) =s2;

  t =s1*conj(s2);
  phi_error =atan(imag(t)/real(t))/(2*pi);
  phi_error_reg(nn) = phi_error;
  sum1 =Kp*phi_error + phi_error*Ki+reg3;

  reg1_reg(nn) =reg1;
  reg2_reg(nn) = reg2;
  reg1 =phi1;

  reg2=reg2+sum1;
  reg3 =reg3+phi_error*Ki;
  phi2_reg(nn) =reg2;
end

figure(1)
plot(phi1_reg);
hold on
plot(phi2_reg,'r');
hold off;
grid on;
title('phase plot');
xlabel('Samples');
ylabel('Phase');

figure(2)

plot(phi_error_reg);
title('phase Error of phase detector');
grid on;
xlabel('samples(n)');
ylabel('Phase error(degrees)');


figure(3)
plot(real(s1_reg));
hold on;
plot(real(s2_reg),'r');
hold off;
grid on;
title('Input signal & Output signal of VCO');
xlabel('Samples');
ylabel('Amplitude');
axis([0 n_data -1.1 1.1]);