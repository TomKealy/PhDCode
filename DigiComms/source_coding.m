function coded_seq=source_coding(src,symbols,codewords)
% Encode a data sequence src based on the given (symbols,codewords).
no_of_symbols=length(symbols);  coded_seq=[];
if length(codewords)<no_of_symbols
  error('The number of codewords must equal that of symbols');
end
for n=1:length(src)
   found=0;
   for i=1:no_of_symbols
      if src(n)==symbols(i), tmp=codewords{i}; found=1; break; end
   end
   if found==0, tmp='?'; end
   coded_seq=[coded_seq tmp];
end