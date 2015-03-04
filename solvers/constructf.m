clear all;
theta = exp(-(1i*2*pi)/(20));
j = -10:10;
j=fliplr(j)
y = zeros(21,21);
for i=1:21
    for k=1:21
        y(i,k) = theta.^(j(k)*i);
    end
end
