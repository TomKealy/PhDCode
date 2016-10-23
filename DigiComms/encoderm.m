function x = encoderm(m,G,map,puncture)
% Copyright 1998, Yufei Wu, MPRG lab, Virginia Tech. for academic use 
% map: Interleaver mapping
% If puncture=0(unpunctured), it operates with a code rate of 1/3. 
% If puncture>0(punctured), it operates with a code rate of 1/2. 
% Multiplexer chooses odd/even-numbered parity bits from RSC1/RSC2.
[N,L] = size(G); % Number of output pits, Constraint length
M = L-1; % Dimension of the state
lm = length(m); % Length of the information message block
lu = lm + M; % Length of the input sequence
% 1st RSC coder output
x1 = rsc_encode(G,m,1);
% interleave input to second encoder
mi = x1(1,map);  x2 = rsc_encode(G,mi,0);
% parallel to serial multiplex to get the output vector
x = [];
if puncture==0   % unpunctured, rate = 1/3;
  for i=1:lu
     x = [x x1(1,i) x1(2,i) x2(2,i)]; 
  end
 else            % punctured into rate 1/2
  for i=1:lu
     if rem(i,2), x = [x x1(1,i) x1(2,i)]; % odd parity bits from RSC1
      else  x = [x x1(1,i) x2(2,i)];       % even parity bits from RSC2    
     end 
  end  
end
x = 2*x - 1; % into bipolar format (+1/-1)
