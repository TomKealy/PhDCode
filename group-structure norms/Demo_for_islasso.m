Seed = 0;
rand('twister',Seed)
randn('state',Seed)

% 1000 data points of dimension 100, on a 10 x 10 grid
n = 1000;
p = 100;

J = [27:29 37:39 47:49]; %p=100


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


lambda = 1;

clear params 

% In order to use ISlasso, we specifiy that the groups are for the
% intersected setting
params.intersection = true;

% Formulation with rectangular groups only
[G,D,IDX] = get_groups( 10, 10, params);

w1 = islasso_SDPT3( XtX_over_n, XtY_over_n, G, D, IDX, lambda, [] );

reshape( w1, 10, 10 )

% Formulation with rectangular groups plus Pi/4 groups

params.slope =[0,1];
[G,D,IDX]    = get_groups( 10, 10, params);


w2 = islasso_SDPT3( XtX_over_n, XtY_over_n, G, D, IDX, lambda, [] );

reshape( w2, 10, 10 )