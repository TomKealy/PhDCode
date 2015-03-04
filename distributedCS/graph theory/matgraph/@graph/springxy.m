function e=springxy(g)
% springxy(g) --- find a spring embedding of g
% 
% This routine is very slow. distxy(g) does a good job and is faster.
% 
% REQUIRE THE OPTIMIZATION TOOLBOX

tic;
n = nv(g);
con = isconnected(g);

if (hasxy(g))
    xy0 = getxy(g);
else
    xy0 = 5*randn(n,2);
end

% opts = optimset('TolX',0.1, 'MaxIter', 50*n,'Display', 'off');
opts = optimset('MaxIter', 10*n,'Display', 'final');

%[xy,e]= fminsearch(@spring_energy, xy0, opts ,g);
%[xy,e]= fminunc(@spring_energy, xy0, opts ,g);

[xy,e] = lsqnonlin(@spring_energy_vec, xy0, [], [], opts);
%[xy,e]= fminunc(@spring_energy, xy0, opts ,g);

embed(g,xy);
toc;


function e = spring_energy(xy,g)

vec = spring_energy_vec(xy);
e = sum(vec.^2);

end

function evec = spring_energy_vec(xy)
%n = nv(g);
evec = zeros(n,n);

for u=1:n-1
    for v=u+1:n
        d = norm(xy(u,:)-xy(v,:));
        if has(g,u,v)
            evec(u,v) = d;
        end
        if ~con
            d = min(d,n/2);
        end
        evec(v,u) = 1/sqrt(d);
    end
end


end


end
