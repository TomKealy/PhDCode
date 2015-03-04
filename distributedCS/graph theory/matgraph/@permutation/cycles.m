function c = cycles(p)
% cycles(p) returns a cell array containing the cycle structure of p

c = {}; 

n = length(p);

if n == 0
    return
end

used = zeros(1,n);
count = 0;
while (sum(used)<n)
    count = count + 1;
    avail = find(used == 0);
    i = avail(1);
    current = i;
    used(i) = 1;
    while(true)
        j = p.array(i);
        if used(j)
            break
        else
            current = [current,j];
            i = j;
            used(i)=1;
        end
    end
    c{count} = current;
end
