function M=combis(N,i)
% Creates an error pattern matrix each row of which is an N-dimensional 
%  vector having ones (representing bit errors) not more than i
M=[]; m1=0;
for n=1:i
   ind= combnk([1:N],n); % nchoosek([1:N],n); 
   for m=1:size(ind,1),  m1=m1+1; M(m1,ind(m,:))=1;  end
end