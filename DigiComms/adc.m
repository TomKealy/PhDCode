function  [d,dcode]=adc(a,b,c,code_table)
% Analog-to-Digital Conversion
% Input : analog signal a, boundary vector b, centroid vector c,
%           and code_table
% Output: quantized samples d and the corresponding code dcode 
N=length(c);
if nargin<4, code_table=[0:N-1]'; end
Na=length(a); % dcode=zeros(Na,size(code_table,2));
for n=1:Na
   I=find(a(n)<b(2:N));
   if ~isempty(I), d(n)=c(I(1)); dcode(n,:)=code_table(I(1),:);
    else           d(n)=c(N); dcode(n,:)=code_table(N,:);
   end    
end
