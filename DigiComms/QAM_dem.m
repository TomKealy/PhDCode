function bitseq=QAM_dem(qamseq,b,bpsym)
%BPSK demodulation
if b==1, bitseq=(qamseq>=0); return; end
%2^b-QAM demodulation
N=length(qamseq);
b1=ceil(b/2); b2=b-b1;
g_code1=2*gray_code(b1)-2^b1+1; g_code2=2*gray_code(b2)-2^b2+1;
tmp1=sum([1:2:2^b1-1].^2)*2^b2; 
tmp2=sum([1:2:2^b2-1].^2)*2^b1;
Kmod=sqrt((tmp1+tmp2)/2/(2^b/4));    % Normalization factor
g_code1=g_code1/Kmod; g_code2=g_code2/Kmod; 
bitseq=[];
for i=1:N
   [emin1,i1]=min(abs(real(qamseq(i))-g_code1));
   [emin2,i2]=min(abs(imag(qamseq(i))-g_code2));
   bitseq=[bitseq deci2bin1(i1-1,b1) deci2bin1(i2-1,b2)];
end
if (nargin>2)
  N = length(bitseq)/bpsym; bitmatrix = reshape(bitseq,bpsym,N).';
  for i=1:N, intseq(i)=bin2deci(bitmatrix(i,:));  end
  bitseq = intseq;
end
