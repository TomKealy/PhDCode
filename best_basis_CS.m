function [xhat ] = best_basis_CS(y, A, mu)

max_iter = 1000;
f=0;
mubar = 1/mu;

fast_sthresh = @(x,th) sign(x).*max(abs(x) - th,0);


for ii = 1:max_iter
   f = f + mubar*A'*(y-A*f);
   
   %minimise the join Largrangian
   T = wpdec(f, 4, 'sym4');
   t = besttree(T);
   
   
   %Denoise by soft thresholding
   f = fast_sthresh(f, )
   
end
xhat = f;
end