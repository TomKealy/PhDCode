M = 300;
K = 100;

edges = [50, 120, 170, 192, 220, 224, 256, 300] ;
levels = [200,  0 , 400, 0, 0, 200, 0, 0];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
g = levels(cumsum(idxs)+1);

igma_squared = [1, 10, 100, 1000];
noise = normrnd(0, 1, [1, M]);
runs = 1;% size(sigma_squared,2);
max_not_in_positions = 0;
ge = g + noise;

% A template is given
temp = [zeros(100,1) ; 100*ones(100,1) ; zeros(100,1)];

% Create a matched filter based on the template
b = flipud(temp(:));

% For testing the matched filter, create a random signal which
% contains a match for the template at some time index
x = [randn(200,1); ge' ; randn(300,1)];
n = 1:length(x);

% Process the signal with the matched filter
y = filter(b,1,x);

% Set a detection threshold (exmaple used is 90% of template)
thresh = 0.9

% Compute normalizing factor
u = temp.'*temp;

% Find matches
matches = n(y>thresh*u);

% Plot the results
plot(n,y,'b', n(matches), y(matches), 'ro');

% Print the results to the console
display(matches);

A = normrnd(0, 1/M, [K, M]);
    
B = diag(diag(A'*A));
    
C = A'*A - B;

y = A*g';

t = A*temp;

b = flipud(t(:));

z = filter(b, 1, y);

% Set a detection threshold (exmaple used is 90% of template)
thresh = 0.9

% Compute normalizing factor
u = t.'*t;

n =1:100;

% Find matches
matches = n(z>thresh*u);

% Plot the results
figure
plot(n,z,'b', n(matches), z(matches), 'ro');

% Print the results to the console
display(matches);


    
