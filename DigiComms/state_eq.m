function [nxb,yb]=state_eq(xb,u,G) 
% To be used as a subroutine for conv_encoder()
K=length(u); LK=size(G,2); L1K=LK-K;
if isempty(xb), xb=zeros(1,L1K); 
 else 
  N=length(xb); %(L-1)K
  if L1K~=N, error('Incompatible Dimension in state_eq()'); end
end
A=[zeros(K,L1K); eye(L1K-K) zeros(L1K-K,K)]; 
B=[eye(K); zeros(L1K-K,K)];
C=G(:,K+1:end); D=G(:,1:K);
nxb=rem(A*xb'+B*u',2)';
yb=rem(C*xb'+D*u',2)';
