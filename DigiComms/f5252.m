function y = f5252(x,SNRb,b)   		
y=Q(-sqrt(2)*x-sqrt(b*SNRb)).^(2^b-1); % Eq.(5.2.52) except e^(-x^2)
