function [decoded_seq,dictionary]=LZW_decoding(c,dictionary)
% Lempel-Ziv-Welch decoding procedure 
% c: LZW-coded sequence in a string
if nargin<2, dictionary={'0','1'}; end
codeset='0123456789abcdefghijklmnopqrstuvwxyz';
lc=length(c);
decoded_seq=[]; decode=[]; 
for n=1:lc
   [decode,dictionary]=LZW_decode(c(n),dictionary,decode,codeset);
   decoded_seq=[decoded_seq decode];
end
