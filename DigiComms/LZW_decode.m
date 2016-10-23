function [decode,dictionary]=LZW_decode(cn,dictionary,pdp,codeset)
% Lempel-Ziv-Welch decoding algorithm 
% cn        = Recent LAW code in character
% dictionary= Set of phrases as a cell like {'0','1','01'}
% pdp       = Previous decoded phrase 
decode_no=findstr(codeset,cn);
if decode_no<=length(dictionary) % Already in the dictionary
  decode=dictionary{decode_no};
  candidate=[pdp decode(1)]; 
 else % Not in the existent dictionary
  candidate=[pdp pdp(1)];
  % Candidate for a new phrase to add into the dictionary
  decode=candidate;
end
if isempty(strmatch(candidate,dictionary,'exact'))
  dictionary{length(dictionary)+1}=candidate; % Add a new phrase
end