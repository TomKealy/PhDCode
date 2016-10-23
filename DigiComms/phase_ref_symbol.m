function PRS=phase_ref_symbol()
% Nfft=Nsd+Nvc=1536+512=2048; 
Nfft=1536+512; %Nsd(# of data subcarriers)+Nvc(# of virtial carrier)
h= ...
[0 2 0 0  0 0 1 1  2 0 0 0  2 2 1 1  0 2 0 0  0 0 1 1  2 0 0 0  2 2 1 1;
 0 3 2 3  0 1 3 0  2 1 2 3  2 3 3 0  0 3 2 3  0 1 3 0  2 1 2 3  2 3 3 0;
 0 0 0 2  0 2 1 3  2 2 0 2  2 0 1 3  0 0 0 2  0 2 1 3  2 2 0 2  2 0 1 3;
 0 1 2 1  0 3 3 2  2 3 2 1  2 1 3 2  0 1 2 1  0 3 3 2  2 3 2 1  2 1 3 2;
0 2 0 0  0 0 1 1  2 0 0 0  2 2 1 1  0 2 0 0  0 0 1 1  2 0 0 0  2 2 1 1];
n=[1 2 0 1  3 2 2 3  2 1 2 3  1 2 3 3  2 2 2 1  1 3 1 2 ...
   3 1 1 1  2 2 1 0  2 2 3 3  0 2 1 3  3 3 3 0  3 0 1 1];
jpi2 = j*pi/2;
for p=1:12 %12*4*32=1536
  for i=1:4
    if p<=6, i1=i; else i1=6-i; end
   for k=1:32 
       temp_seq((p-1)*128+(i-1)*32+k)=exp(jpi2*(h(i1,k)+n((p-1)*4+i))); 
    end
  end
end
PRS = [zeros(1,256) temp_seq(1:768) 0 temp_seq(769:end) zeros(1,255)];