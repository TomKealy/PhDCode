function X = pm_haar(x,s)
%
% MATLAB function for pm Haar 1-D DWT algorithm. 
% This function receives N (a power of two) real values in x
% and the number of stages (scales) of decomposition required in s, 
% computes the Haar DWT, and returns the N DWT coefficients in X.
%
% For the theory behind this algorithm and example input/output,
% please refer the paper:
% Fundamentals of the discrete Haar wavelet transform
% Duraisamy Sundararajan
% dsprelated.com, 2011
% articles/paper section
%
N = length(x);  % length of the input vector 
b = N;sq2 = sqrt(2);fac = 1;
  for ns =1:s % outer loop stepping over stages
    p = 1;q = 1;
    while(p < b) % inner loop stepping over each butterfly in a stage
      r = p + 1;
      xmr(q) = x(p) - x(r); xpr(q) = x(p) + x(r);
      p = p + 2; q = q + 1;
    end
    fac = fac / sq2;
    x(1:b) = [xpr(1:b/2) xmr(1:b/2) * fac];
    b = b / 2;
  end
x(1:b) = fac * x(1:b);
X=x;