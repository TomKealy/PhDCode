M = zeros(10,10)

for i=1:10
    for j=1:10
        M(i,j) = min(i,j);
    end
end
