function [h,L,H]=Huffman_code(p,opt)
%Huffman code generator gives a Huffman code matrix h, 
%  average codeword length L & entropy H
% for a source with probability vector p given as argin. 
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
zero_one=['0';'1']; 
if nargin>1&&opt>0, zero_one=['1';'0']; end
if abs(sum(p)-1)>1e-6
  fprintf('\n The probabilities in huffman_code() does not add up to 1!');
  p=p/sum(p);
end  
M=length(p);  N=M-1; p=p(:); % make p a column vector
h{1}=zero_one(1); h{2}=zero_one(2);
h={zero_one(1),zero_one(2)};
if M>2
  pp(:,1)=p;
  for n=1:N
     [pp(1:M-n+1,n),o(1:M-n+1,n)]=sort(pp(1:M-n+1,n),1,'descend'); % in descending order
     if n==1, ord0=o; end  % Original descending order
     if M-n>1
       pp(1:M-n,n+1)=[pp(1:M-1-n,n); sum(pp(M-n:M-n+1,n))];
     end
  end
  for n=N:-1:2
     tmp=N-n+2; oi=o(1:tmp,n);
     for i=1:tmp, h1{oi(i)}=h{i}; end
     h=h1; h{tmp+1}=h{tmp};
     h{tmp}=[h{tmp} zero_one(1)]; 
     h{tmp+1}=[h{tmp+1} zero_one(2)];
  end
  for i=1:length(ord0), h1{ord0(i)}=h{i}; end
  h=h1;
end
L=0; 
for n=1:M
  L=L+p(n)*length(h{n});% Average codeword length
end
H=-sum(p.*log2(p)); % Entropy by Eq.(9.1.4)