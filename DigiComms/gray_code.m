function [g_code,b_code]=gray_code(b)
N=2^b; g_code=0:N-1; 
if b>1, g_code=gray_code0(g_code); end
b_code=deci2bin(g_code);