function L_A = sova(Ly,G,Lu,ind_dec) 
% Copyright: Yufei Wu, 1998, MPRG lab, Virginia Tech for academic use
% This implements Soft Output Viterbi Algorithm in trace back mode 
% Input: Ly     : Scaled received bits Ly=0.5*L_c*y=(2*a*rate*Eb/N0)*y
%        G      : Code generator for the RSC code in binary matrix form
%        Lu     : Extrinsic information from the previous decoder.
%        ind_dec: Index of decoder=1/2
%                  (assumed to be terminated in all-zero state/open)
% Output: L_A  : Log-Likelihood Ratio (soft-value) of 
%                   estimated message input bit u(k) at each stage, 
%                 ln (P(u(k)=1|y)/P(u(k)=-1|y))
lu = length(Ly)/2; % Number of y=[ys yp] in Ly
lu1 = lu+1; Infty = 1e2;
[N,L] = size(G); Ns = 2^(L-1); % Number of states
delta = 30;  % SOVA window size
%Make decision after 'delta' delay. Tracing back from (k+delta) to k,
% decide bit k when received bits for bit (k+delta) are processed.  
% Set up the trellis defined by G.
[nout,ns,pout,ps] = trellis(G);
% Initialize the path metrics to -Infty
Mk(1:Ns,1:lu1)=-Infty; Mk(1,1)=0; % Only initial all-0 state possible
% Trace forward to compute all the path metrics
for k=1:lu
   Lyk = Ly(k*2-[1 0]); k1=k+1;
   for s=1:Ns  % Eq.(9.4.44), Eq.(9.4.45)
      Mk0 = Lyk*pout(s,1:2).' -Lu(k)/2 +Mk(ps(s,1),k);
      Mk1 = Lyk*pout(s,3:4).' +Lu(k)/2 +Mk(ps(s,2),k);
      if Mk0>Mk1, Mk(s,k1)=Mk0; DM(s,k1)=Mk0-Mk1; pinput(s,k1)=0;
       else       Mk(s,k1)=Mk1; DM(s,k1)=Mk1-Mk0; pinput(s,k1)=1;
      end
   end
end
% Trace back from all-zero state or the most likely state for D1/D2
% to get input estimates uhat(k), and the most likely path (state) shat
if ind_dec==1, shat(lu1)=1;  else [Max,shat(lu1)]=max(Mk(:,lu1));  end
for k=lu:-1:1
   uhat(k)=pinput(shat(k+1),k+1); shat(k)=ps(shat(k+1),uhat(k)+1);
end
% As the soft-output, find the minimum DM over a competing path 
%  with different information bit estimate.
for k=1:lu
   LLR = min(Infty,DM(shat(k+1),k+1));
   for i=1:delta
      if k+i<lu1
        u_=1-uhat(k+i); % the information bit
        tmp_state = ps(shat(k+i+1),u_+1);
        for j=i-1:-1:0
           pu=pinput(tmp_state,k+j+1); tmp_state=ps(temp_state,pu+1);
        end
        if pu~=uhat(k), LLR = min(LLR,DM(shat(k+i+1),k+i+1)); end
      end
   end
   L_A(k) = (2*uhat(k)-1)*LLR; % Eq.(9.4.46)
end