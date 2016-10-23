function qamseq=QAM(bitseq,b)
bpsym = nextpow2(max(bitseq)); % no of bits per symbol
if bpsym>0, bitseq = deci2bin(bitseq,bpsym);  end
if b==1, qamseq=bitseq*2-1; return; end  % BPSK modulation
% 2^b-QAM modulation
N0=length(bitseq);  N=ceil(N0/b);
bitseq=bitseq(:).'; bitseq=[bitseq zeros(1,N*b-N0)];
b1=ceil(b/2); b2=b-b1; b21=b^2; b12=2^b1; b22=2^b2; 
g_code1=2*gray_code(b1)-b12+1; g_code2=2*gray_code(b2)-b22+1;
tmp1=sum([1:2:2^b1-1].^2)*b21; tmp2=sum([1:2:b22-1].^2)*b12;
M=2^b; Kmod=sqrt(2*(M-1)/3); 
%Kmod=sqrt((tmp1+tmp2)/2/(2^b/4)) % Normalization factor
qamseq=[];
for i=0:N-1
   bi=b*i; i_real=bin2deci(bitseq(bi+[1:b1]))+1;
   i_imag=bin2deci(bitseq(bi+[b1+1:b]))+1;
   qamseq=[qamseq (g_code1(i_real)+j*g_code2(i_imag))/Kmod];
end