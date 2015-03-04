function [delays,distances] = generate_delays(xy)
%Generates an array of delays from (0,0) to (x,y) by calculating the
%euclidiean distance and dividing by c
[num_nodes,dim] = size(xy);
s = zeros(1,num_nodes);
for i=1:num_nodes
    s(i) = norm(xy(i,:),2);
end
c = 3e8;
delays = s/c;
distances = s;
end

