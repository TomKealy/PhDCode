function AI=inv_GF2(A)
[M,N]=size(A);  I=eye(M);
for i=1:M, AI(:,i)=GFlineq(A,I(:,i)); end