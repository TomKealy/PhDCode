function y = proj_group_20110627(z,w,nonneg,groups,g)
%
%  Group-wise or row-wise projection operator
%

if nonneg  % nonnegativity
    zn = min(z,0);
    z = max(z,0);
end

if nargin == 5  % group-wise
    t = min(1,w./sqrt(g*z.^2));
    y = t(groups).*z;
end

if nargin == 3  % row-wise
    t = min(1,w./sqrt(sum(z.^2,2)));
    y = bsxfun(@times,z,t);
end

if nonneg  % nonnegativity
    y = max(y,0) + zn;
end

end