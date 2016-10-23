function [y,ch_buf] = channel(x,h,ch_buf)
L_x = length(x); L_h = length(h); h=h(:); 
if nargin<3, ch_buf = zeros(1,L_h); 
 else L_ch_buf = length(ch_buf);
      if L_ch_buf<L_h,  ch_buf = [ch_buf zeros(1,L_h-L_h)]; 
       else  ch_buf = ch_buf(1:L_h);
      end    
end
for n=1:L_x,  ch_buf=[x(n) ch_buf(1:end-1)];  y(n)=ch_buf*h;  end