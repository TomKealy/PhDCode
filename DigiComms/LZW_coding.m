function coded_seq=LZW_coding(src,dictionary)
% Lempel-Ziv-Welch encoding procedure 
% src: Source code in a binary number string
ls=length(source);
src_str=[]; % Initialize message string
if ~isstr(src) % If the source code consists of numbers, not string
  for i=1:ls, src_str=[src num2str(src(i))]; end
  src=src_str;
end
w=[]; coded_seq=[];
codeset='0123456789abcdefghijklmnopqrstuvwxyz';
for n=1:ls
   [code,dictionary,w]=LZW_code(src(n),dictionary,w,codeset);
   coded_seq=[coded_seq code];
end
code=codeset(strmatch(w,dictionary,'exact'));
coded_seq=[coded_seq code];