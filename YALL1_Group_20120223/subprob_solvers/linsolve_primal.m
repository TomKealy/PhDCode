function [x,Ax] = linsolve_primal(A,At,b,z,beta,lambda1,lambda2,opts,invIpAAt,x,Ax)
%
%  Solve the linear system:
%
%    (beta1*I+beta2*A'*A)x = beta1*z-lambda1+A'*(beta2*b+lambda2),
%
%  either exactly or approximately.
%
%  Ax = A*x is saved for other uses in the main solver
%
%  To save computing time in case of non-orthonormal A, exact computation
%  uses operator invIpAAt (see YALL1_group.m), and inexact computation uses
%  gradient descent and takes advantages of existing x and Ax

if ~opts.nonorth  % orthonormal A
    u = beta(1)*z - lambda1;
    v = beta(2)*b + lambda2;
    Au = A(u);
    x = u/beta(1) - At(beta(2)/beta(1)*Au - v)/(beta(1)+beta(2));
    Ax = (Au + v)/(beta(1)+beta(2));
else  % non-orthonormal A
    if opts.exact  % solve a linear system exactly
        x = beta(1)*z - lambda1 + At(beta(2)*b + lambda2);
        x = 1/beta(1)*(x - beta(2)*At(invIpAAt(A(x))));
        Ax = A(x);
    else           % take a gradient descent step
        grad = beta(1)*(x-z)+lambda1+At(beta(2)*(Ax-b)-lambda2);  % gradient
        Ag = A(grad);
        gtg = sum(grad.^2);
        stp = gtg./(beta(1)*gtg+beta(2)*sum(Ag.^2)); % step length
        x = x-bsxfun(@times,grad,stp);
        Ax = Ax-bsxfun(@times,Ag,stp);
    end
end