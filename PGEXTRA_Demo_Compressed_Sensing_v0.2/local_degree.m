function W_FDLA=local_degree(P,eps_deg)
if ~exist('eps_deg','var')
    eps_deg=1;
end;
L=size(P,1);
W_FDLA(1:L,1:L)=0;
deg=sum(P,2);
for i=1:L
    for j=1:L
        if P(i,j)==1
            W_FDLA(i,j)=1/(max(deg(i),deg(j))+eps_deg);
        end;
    end;
end;
W_FDLA=diag(ones(L,1)-sum(W_FDLA,2))+W_FDLA;
end