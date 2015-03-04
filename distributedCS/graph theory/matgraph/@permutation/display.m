function st = display(p)
% display(p) --- display a permutation in disjoint cycle form. 


% special case if the permutation is empty 
if (length(p)==0)
    disp('()')
    return
end

c = cycles(p);
outstr = '';

for k=1:length(c)
    outstr = [outstr,list_to_cycle(c{k})];
end


if nargout > 0
	st = outstr;
else
	disp(outstr)
end


function str = list_to_cycle(list)
if length(list) == 0
    str = '()'
    return
end

str = '(';
for k=1:length(list)
    str = [str,int2str(list(k))];
    if k < length(list)
        str = [str,','];
    else
        str = [str,')'];
    end
end