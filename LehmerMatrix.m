function [M] = LehmerMatrix( n )

M = zeros(n,n);

for i=1:n
    for j=1:n
        M(i,j) = min(i,j);
    end
end

end

