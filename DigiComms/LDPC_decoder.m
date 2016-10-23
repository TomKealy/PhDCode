function vhat=LDPC_decoder(H,Lcy)
% Performs an LDPC decoding with given parity-check matrix H and Lc*y 
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
[M,N] = size(H); %rate = M/N;
Lcy = Lcy(1:N); 
% Prior Probabilities
pn0=1./(1+exp(Lcy)); pn1=1-pn0; % for u(i)=-1 (m=0)/u(i)=1 (m=1)
% Initialization
dv=max(sum(H)); dc=max(sum(H.')); Nm=zeros(M,dc); Mn=zeros(dv,N);
for m=1:M
   ind = 1; 
   for n=1:N
      if H(m,n)==1
        Nm(m,ind)=n; ind =ind+1; qmn0(m,n)=pn0(n); qmn1(m,n)=pn1(n); 
      end
   end
end
for n=1:N
   tmp = find(H(:,n)); lt=length(tmp); if lt>0, Mn(1:lt,n)=tmp; end
end
for iter=1:40
   for m=1:M % Horizontal step for every check node
      NNm = sum(Nm(m,:)~=0);
      for n0=1:NNm,  n=Nm(m,n0); dqmn(m,n0)=qmn0(m,n)-qmn1(m,n);  end
      for n0=1:NNm
         n = Nm(m,n0); drmn=1; tmp0=1; tmp1=1;
         for k=1:NNm,  if k~=n0, drmn=drmn*dqmn(m,k); end,  end
         rmn0(m,n)=(1+drmn)/2; rmn1(m,n)=1-rmn0(m,n); % Eq.(9.4.56)
      end
   end % End of horizontal step
   for n=1:N % Vertical step for every variable node
      NMn = sum(Mn(:,n)~=0); 
      for m0=1:NMn
         m = Mn(m0,n); pn0_prod_rmn0=pn0(n); pn1_prod_rmn1=pn1(n);
         for k=1:NMn
            if k~=m0
              pn0_prod_rmn0 = pn0_prod_rmn0*rmn0(m,n);
              pn1_prod_rmn1 = pn1_prod_rmn1*rmn1(m,n);
            end
         end 
         alphamn=1/(pn0_prod_rmn0+pn1_prod_rmn1);   
         qmn0(m,n)=alphamn*pn0_prod_rmn0; % Eq.(9.4.57)
         qmn1(m,n)=1-qmn0(m,n);  
      end
      %update pseudo posterior probability
      pn0_prod_all_rmn0 = pn0_prod_rmn0*rmn0(m,n);
      pn1_prod_all_rmn1 = pn1_prod_rmn1*rmn1(m,n);
      alpha_n=1/(pn0_prod_all_rmn0+pn1_prod_all_rmn1);
      qn1(n) = alpha_n*pn1_prod_all_rmn1; % Eq.(9.4.58b)
      vhat(n)= qn1(n)>0.5; % % Tentative decoding by Eq.(9.4.59)   
   end % End of vertical step
   if sum(rem(vhat*H.',2))==0 | (iter>20 & vhat==vh0), break; end
   vh0 = vhat;
end % End of iter loop