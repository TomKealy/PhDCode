function [H,G]=Hamm_gen(m,opt)
% Generates the parity-check/generator matrices for Hamming Code.
% Output: H= mxN parity check matrix (N=2^m-1)
%         G= KxN generator matrix for Hamming Code (K=N-m)
if nargin<2,  opt=0;  end
H=[ones(2,1) eye(2)];  
G=[1 H(:,1).'];
if m<3,  return;  end
for i=3:m
   N=2^i; N2=N/2; K=N-1-i;  
   for j=1:N2-1, H(i,j) = 1;  end
   for j=N2:K,  H(i,j) = 0;  end
   H(1:i-1,N2:K) = H(1:i-1,1:N2-i);
   H(1:i,N-i:N-1) = eye(i);
end
H(:,1:K) = fliplr(H(:,1:K));
if opt>0
  G = [eye(K) H(1:m,1:K).'];
 else
  H = [H(:,N-m:N-1) H(:,1:N-m-1)];
  G = [H(1:i,i+1:N-1).' eye(K)];
end