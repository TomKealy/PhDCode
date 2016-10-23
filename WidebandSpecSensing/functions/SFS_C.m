% find sample pattern with sequential forward search
function [C,KA]=SFS_C(L,p,S,Fs)
% L period ofcoset sampling
% p number of coset samples
% S spectral set: S is subset of {0 to L-1}
% Fs Nyquist frequency
%% SFS pattern
T=1/Fs;
k=1:L;
KA=[];
ka1=[];
fin_KA=[];
C_SFS=[];
for st=0:0
    C_ls=[st];
    C1=0:L-1;
    ka1=[];
    for m=1:p-1
        KA=[];
        C1=setdiff(C1,C_ls);
        for h=1:length(C1)
            C=[C_ls,C1(h)];
            C=sort(C);
            A=1/(L*T)*exp(i*2*pi*C'*(k-1)/L);
            As=A(:,S);
            KA=[KA,cond(As)];
        end
        [min_KA cmin]=min(KA);
        C_ls=[C_ls,C1(cmin)];
        ka1=[ka1,min_KA];
    end
C_ls=sort(C_ls);
C=C_ls;
A=1/(L*T)*exp(i*2*pi*C'*(k-1)/L);
As=A(:,S);
fin_KA=[fin_KA;cond(As)];
KA_SFS=ka1;
C_SFS=[C_SFS;C];
end
C=C_SFS;
KA=fin_KA;
end