function display(p)
% display(p) --- used to display paritions to the console

[m,n] = size(p.array);

outstr = '{ ';

for k=m:-1:1
    part = find(p.array(k,:));
    outstr = [outstr,disp_part(part),' '];
end
outstr = [outstr,'}'];
disp(outstr);

end



function s = disp_part(part)
s = '{';
for k=1:length(part)
    s = [s,int2str(part(k))];
    if (k<length(part))
        s = [s,','];
    end
end
s = [s,'}'];
end
