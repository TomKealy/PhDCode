function [x,Ax,Gx] = linsolve_overlap_primal_20120219(A,At,G,b,z,beta,lambda1,lambda2,x,Ax,Gx)
%
%  Solve the linear system:
%
%    (beta1*G'*G+beta2*A'*A)x = G'*(beta1*z-lambda1)+A'*(beta2*b+lambda2),
%
%  approximately by one-step gradient descent.
%
%  Ax = A*x, Gx = G*x are saved for other uses in the main solver.
%------------------------------------------------------------------------

grad = ((beta(1)*(Gx-z)+lambda1)'*G)'+At(beta(2)*(Ax-b)-lambda2);  % gradient
Ag = A(grad); Gg = G*grad;
stp = sum(grad.^2)/(beta(1)*sum(Gg.^2)+beta(2)*sum(Ag.^2)); % step length
x = x-stp*grad;
Ax = Ax-stp*Ag;
Gx = Gx-stp*Gg;