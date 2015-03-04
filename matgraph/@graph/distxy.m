function e = distxy(g, D) 
% distxy(g) -- give g a distance based embedding
% we attempt to embed g in the plane so that the graph-theoretic distance
% between vertices matches the eucliden distance.
% 
% This may also be called distxy(g,D) where D is a distance matrix. (By
% default D is the standard shortest-path distance, but could be replaced
% by, say, resistance(g).
% 
% REQUIRES THE OPTIMIZATION TOOLBOX

tic;
n = nv(g);

if nargin ==1
	d = dist(g);
else
	d = D;
end

[i,j] = find(d==inf);
ni = length(i);
for k=1:ni
    d(i(k),j(k)) = n/2;
end


if (hasxy(g))
    xy0 = getxy(g);
else
    xy0 = 5*randn(n,2);
end

opts = optimset('MaxIter', 5*n,'Display', 'final');

[xy,e] = lsqnonlin(@dist_discrep, xy0, [], [], opts);
disp(['Embedding score = ', num2str(e)])
embed(g,xy);
toc

function dd = dist_discrep(xy)

%divisor = d + eye(n);
dxy = dist_pair(xy);
divisor = real((d+eye(n)));
% divisor = exp(d/4)-0.8;
dd = d - dxy;
dd = dd./divisor;
dd = dd(:);

end
    
    
%--------------------------------------------------------
function D = dist_pair(xy)
% find Euclidean distances between vertices in xy embedding
n = size(xy,1); % get # of rows
c = sum(xy.*xy,2); % norm^2 of rows
Y = c * ones(1,n);
D = Y + Y' - 2*xy*xy';
D = real(sqrt(D));

end % end of dist_pair

end % end of distxy
