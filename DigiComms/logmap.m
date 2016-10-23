function L_A = logmap(Ly,G,Lu,ind_dec)
% Copyright 1998, Yufei Wu, MPRG lab, Virginia Tech. for academic use
% Log_MAP algorithm using straightforward method to compute branch cost
% Input: Ly     = scaled received bits Ly=0.5*L_c*y=(2*a*rate*Eb/N0)*y
%        G      = code generator for the RSC code a in binary matrix 
%        Lu     = extrinsic information from the previous decoder.
%        ind_dec= index of decoder=1/2 (assumed to be terminated/open)
% Output: L_A   = ln (P(x=1|y)/P(x=-1|y)), i.e., Log-Likelihood Ratio 
%              (soft-value) of estimated message input bit at each level 
lu=length(Ly)/2; Infty=1e2; EPS=1e-50; % Number of input bits, etc
[N,L] = size(G); 
Ns = 2^(L-1);  % Number of states in the trellis
Le1=-log(1+exp(Lu)); Le2=Lu+Le1; % ln(exp((u+1)/2*Lu)/(1+exp(Lu)))
% Set up the trellis
[nout, ns, pout, ps] = trellis(G);
% Initialization of Alpha and Beta
Alpha(1,2:Ns) = -Infty; % Eq.(9.4.40) (the initial all-zero state) 
if ind_dec==1 % for decoder D1 with termination in all-zero state
  Beta(lu+1,2:Ns) = -Infty; % Eq.(9.4.41b) (the final all-zero state)
 else % for decoder D2 without termination
  Beta(lu+1,:) = -log(Ns)*ones(1,Ns);
end
% Compute gamma at every depth level (stage)
for k = 2:lu+1
   Lyk = Ly(k*2-[3 2]);  gam(:,:,k) = -Infty*ones(Ns,Ns);
   for s2 = 1:Ns % Eq.(9.4.42)
      gam(ps(s2,1),s2,k) = Lyk*[-1 pout(s2,2)].' +Le1(k-1);
      gam(ps(s2,2),s2,k) = Lyk*[+1 pout(s2,4)].' +Le2(k-1);
   end
end  
% Compute Alpha in forward recursion
for k = 2:lu
   for s2 = 1:Ns
      alpha = sum(exp(gam(:,s2,k).'+Alpha(k-1,:))); % Eq.(9.4.40)
      if alpha<EPS, Alpha(k,s2)=-Infty; else Alpha(k,s2)=log(alpha); end   
   end
   tempmax(k) = max(Alpha(k,:)); Alpha(k,:) = Alpha(k,:)-tempmax(k);
end     
% Compute Beta in backward recursion
for k = lu:-1:2
   for s1 = 1:Ns
      beta = sum(exp(gam(s1,:,k+1)+Beta(k+1,:))); % Eq.(9.4.41)
      if beta<EPS, Beta(k,s1)=-Infty; else Beta(k,s1)=log(beta); end   
   end
   Beta(k,:) = Beta(k,:) - tempmax(k);
end
% Compute the soft output LLR for the estimated message input
for k = 1:lu
   for s2 = 1:Ns % Eq.(9.4.39)
      temp1(s2)=exp(gam(ps(s2,1),s2,k+1)+Alpha(k,ps(s2,1))+Beta(k+1,s2));
      temp2(s2)=exp(gam(ps(s2,2),s2,k+1)+Alpha(k,ps(s2,2))+Beta(k+1,s2));
   end
   L_A(k) = log(sum(temp2)+EPS) - log(sum(temp1)+EPS); % Eq.(9.4.38)
end