function [decodes,E,epi]=cyclic_decoder(code_seq,N,K,g,E,epi)
% Cyclic (N,K) decoding of received code_seq with generator polynml g
% E:   Error Pattern matrix or syndromes
% epi: error pattern index vector
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
if nargin<6 
  nceb=ceil((N-K)/log2(N+1)); % Number of correctable error bits
  E = combis(N,nceb);  % All error patterns
  for i=1:size(E,1)
     syndrome=cyclic_decoder0(E(i,:),N,K,g); 
synd_decimal=bin2deci(syndrome);
     epi(synd_decimal)=i; % Error pattern indices
  end
end
if (size(code_seq,2)==1)  code_seq=code_seq.';  end
Lcode= length(code_seq);  Ncode= ceil(Lcode/N);
Code_seq= [code_seq(:); zeros(Ncode*N-Lcode,1)];
Code_seq= reshape(Code_seq,N,Ncode).';
decodes=[]; syndromes=[];
for n=1:Ncode
  code= Code_seq(n,:);
  syndrome= cyclic_decoder0(code,N,K,g);
  si= bin2deci(syndrome);  % Syndrome index
  if 0<si&si<=length(epi)  % Syndrome index to error pattern index
    m=epi(si);  if m>0, code=rem(code+E(m,:),2); end  % Eq.(9.4.28)
  end
  decodes=[decodes code(N-K+1:N)]; syndromes=[syndromes syndrome];
end
if nargout==2, E=syndromes; end