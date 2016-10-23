Seed = 0;
rand('twister',Seed)
randn('state',Seed)

% 1000 data points of dimension 100, on a 10 x 10 grid
n = 1000;
p = 100;

% We define a nonzero pattern
J = [27:29 37:39 47:49];


w    = zeros(p,1);
w(J) = randn(length(J),1);    

% Plot the generating pattern
reshape( w, 10, 10 )


X     = randn(n,p);
noise = randn(n,1);
X_w   = X*w;

SNR = 3;

sigma = norm(X_w)/(norm(noise)*SNR);

Y  = X_w+sigma*noise;

XtX_over_n = X'*X/n;
YtY_over_n = Y'*Y/n;
XtY_over_n = X'*Y/n;


[G,D,IDX] = get_groups( 10, 10, []);

lambda = 0.1;


% Formulation with simple penalization
w1 = activesetalgorithm( X, XtY_over_n, G, D, IDX, lambda, []);
w2 = activesetalgorithm_rectangle( X, XtY_over_n, G, D, IDX, lambda, 10, 10, []);


norm( w1 - w2 )


