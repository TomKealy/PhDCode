function punctured = puncture(un_pun, punc_vector)
% Puncture the bits of un_pun as specified by with punc_vector
L_pun_vector=length(punc_vector);
if L_pun_vector==0, punctured=un_pun; return; end
L_un_pun=length(un_pun); 
L_pass=sum(punc_vector~=0);
L_pun=ceil(L_un_pun/L_pun_vector)*L_pass; 
punctured=[];
for i=1:L_un_pun
   pun_index=mod(i,L_pun_vector);
   if pun_index==0, pun_index=L_pun_vector; end
   if punc_vector(pun_index)~=0, punctured=[punctured un_pun(i)]; end
end
punctured=[punctured zeros(1,L_pun-length(punctured))];
if size(punc_vector,1)>size(punc_vector,2)
  punctured=punctured.';
end