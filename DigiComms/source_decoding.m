function decoded_seq=source_decoding(coded_seq,codewords,symbols)
% Decode a coded_seq  based on the given (codewords,symbols).
M=length(codewords);  decoded_seq=[];
while ~isempty(coded_seq)
  lcs= length(coded_seq); found=0;
  for m=1:M
    codeword= codewords{m};
    lc= length(codeword);
    if lcs>=lc&codeword==coded_seq(1:lc)
      symbol=symbols(m); found=1;  break;
    end
    if found==0, symbol='?'; end
  end 
  decoded_seq=[decoded_seq symbol]; 
  coded_seq=coded_seq(lc+1:end); 
end