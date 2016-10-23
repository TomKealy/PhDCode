function [code,dictionary,w]=LZW_code(dn,dictionary,w,codeset)
% LZW coding algorithm 
% dn        = Recent source code in character
% dictionary= Set of phrases as a cell like {'0','1','01'}
% w         = Stem of phrase 
wdn= [w dn]; % Possible candidate for a new phrase
if isempty(strmatch(wdn,dictionary,'exact')) %not in the dictionary
  code=codeset(strmatch(w,dictionary,'exact'));
  dictionary{length(dictionary)+1}=wdn; %add a new phrase
  w=dn;
 else % Already in the dictionary
  code=[];
  w=wdn; 
end