function B=mulmd(A,d)
%B=mulmd(A,d)
%
% B=A*diag(d)

[M,N]=size(A);
n=length(d);
if n ~= N, error('Wrong matrix sizes');end

B=A;
for i=1:N
   B(:,i)=B(:,i)*d(i);
end

