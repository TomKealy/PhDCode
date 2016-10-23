function c=PNG(gm,Kc,x)
%generates an m-sequence from a feedback-shift register if Kc=0
% or from a high-speed feedback shift generator if KC>0
% gm: a generator tap coefficient vector 
%      or a set of connected tap position numbers
if sum(gm>1)>0, g(gm)=1;  % say, if gm=[1 2 4 5]
 else g=gm;                  % if gm=[1 1 0 1 1]
end
m=length(g); % Number of feedback shift-registers
if nargin==2
  if length(Kc)>1, x=Kc; Kc=0; % If x0 is given as 2nd input argument
   else x=[zeros(1,m-1) 1]; % Initial state of the shift register
  end     
elseif nargin<2, Kc=0; x=[zeros(1,m-1) 1];
end   
x=x(end:-1:1); N=2^m-1; 
if Kc==0 % Feedback shift register
  g=fliplr(g); for i=1:N, c(i)=x(1); x=[x(2:m) rem(g*x',2)]; end
 elseif Kc<0
  g=fliplr(g); c=x;
  for i=m+1:N, x=[x(2:m) rem(g*x',2)]; c(i)=x(m); end
 else % High-speed feedback shift generator
  c=x;
  for i=m+1:N, c(i)=x(1); x=[rem(x(2:m)+g(1:m-1)*x(1),2) x(1)];  end
end