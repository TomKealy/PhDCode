function I=Gauss_Hermite(f,N,varargin)
% Refer to 'Applied Numerical Methods Using MATLAB' by Won Yang et. al
[t,w]=Gausshp(N);
ft=feval(f,t,varargin{:}); 
I= w*ft';
