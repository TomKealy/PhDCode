function [yden,lambda,MSE] = oracle(denoisefun,ytrue,ynoisy,lambdamin,lambdamax)

relerr = 1/100;

%MSEin  = norm(ynoisy-ytrue).^2;
lambda = fmin(@(l) norm(denoisefun(ynoisy,l)-ytrue),lambdamin,lambdamax,relerr);
yden   = denoisefun(ynoisy,lambda);
MSE    = norm(yden-ytrue).^2;

%if MSE > MSEin
%    disp('something is TERRIBLY WRONG!\n');
%end;

function lm = fmin(objfun,lmin,lmax,relerr)

l = exp(linspace(log(lmin),log(lmax),4));
o = [objfun(l(1)) objfun(l(2)) objfun(l(3)) objfun(l(4))];

%fprintf('oracle: start');

[om im] = min(o);
while true
    %fprintf('.');
    switch im
    case 1
        lmin = l(1);
        lmax = l(2);
        l = exp(linspace(log(lmin),log(lmax),4));
        o = [o(1) objfun(l(2)) objfun(l(3)) o(2)]; 
    case 2
        lmin = l(1);
        lmax = l(3);
        l = exp(linspace(log(lmin),log(lmax),4));
        o = [o(1) objfun(l(2)) objfun(l(3)) o(3)]; 
%       l = [lmin (lmin+l(2))/2 l(2) lmax];
%       o = [o(1) objfun(l(2)) o(2) o(3)]; 
%       l = [exp(linspace(log(lmin),log(l(2)),3)) lmax];
%       o = [o(1) objfun(l(2)) o(2) o(3)]; 
    case 3
        lmin = l(2);
        lmax = l(4);
        l = exp(linspace(log(lmin),log(lmax),4));
        o = [o(2) objfun(l(2)) objfun(l(3)) o(4)]; 
%       l = [lmin (lmin+l(3))/2 l(3) lmax];
%       o = [o(2) objfun(l(2)) o(3) o(4)]; 
%       l = [exp(linspace(log(lmin),log(l(3)),3)) lmax];
%       o = [o(2) objfun(l(2)) o(3) o(4)]; 
    case 4
        lmin = l(3);
        lmax = l(4);
        l = exp(linspace(log(lmin),log(lmax),4));
        o = [o(3) objfun(l(2)) objfun(l(3)) o(4)]; 
    end;

    [om im] = min(o);
    if log(abs(lmax)) - log(abs(lmin)) < relerr
        break;
    end;
end;
%fprintf('end\n');

lm = l(im);
