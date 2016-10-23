function [y,Aty] = linsolve_dual_20110620(A,At,b,x,z,beta,opts,invAAt,y,Aty)
%
%  Solve the linear system:
%
%              (beta*A'*A)y = b-A*x+beta*A*z,
%
%  either exactly or approximately.
%
%  To save computing time in case of non-orthonormal A, exact computation
%  uses operator invAAt (see YALL1_group.m), and inexact computation applies
%  gradient descent and takes advantages of existing y and Aty=A'*y

if ~opts.nonorth  % orthonormal A
    y = b/beta + A((z-x/beta));
    Aty = At(y);
else  % non-orthonormal A
    if opts.exact,  % solve a linear system
        y = invAAt((b/beta + A((z-x/beta))));
        Aty = At(y);
    else            % take a gradient descent step
        grad = A(x/beta+Aty-z)-b/beta;
        Atg = At(grad);
        stp = sum(grad.^2)/sum(Atg.^2);  % step length
        y = y-bsxfun(@times,grad,stp);
        Aty = Aty-bsxfun(@times,Atg,stp);
    end
end
