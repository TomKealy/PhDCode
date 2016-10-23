function y = demultiplex(r,map,puncture)
%Copyright 1998, Yufei Wu, MPRG lab, Virginia Tech. for academic use
% map: Interleaver mapping 
Nb = 3-puncture; lu = length(r)/Nb;
if puncture==0   % unpunctured
  for i=1:lu, y(:,2*i) = r(3*i-[1 0]).'; end
else             % punctured
  for i=1:lu
     i2 = i*2;
     if rem(i,2)>0, y(:,i2)=[r(i2); 0]; else y(:,i2)=[0; r(i2)];  end      
  end
end
sys_bit_seq = r(1,1:Nb:end); % the systematic bits for both decoders
y(:,1:2:lu*2) = [sys_bit_seq; sys_bit_seq(map)];