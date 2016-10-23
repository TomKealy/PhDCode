function depunctured = depuncture(pun,punc_vector)
% Insert zeros into the punctured position of pun as specified by punc_vector
L_pun_vector=length(punc_vector);
if L_pun_vector==0, depunctured=pun; return; end
L_pun=length(pun);
L_zero=sum(punc_vector==0);
L_depun=L_pun+floor(L_pun/(L_pun_vector-L_zero))*L_zero;
depunctured=zeros(1,L_depun);
ip=1;  id=1;
while ip<=L_pun
   ipv=mod(id,L_pun_vector);
   if ipv==0, ipv=L_pun_vector;  end
   if punc_vector(ipv)~=0
     depunctured(id)=pun(ip); ip=ip+1; 
   end
   id=id+1;
end
if size(punc_vector,1)>size(punc_vector,2)
  depunctured=depunctured.';  
end