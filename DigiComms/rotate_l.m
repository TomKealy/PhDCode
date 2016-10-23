function xl=rotate_l(x,M)
N=size(x,2); M=mod(M,N); xl=[x(:,M+1:end) x(:,1:M)];