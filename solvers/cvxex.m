clear
% size of signal
n = 512;
x0 = zeros(n,1);
% function is f(x)=2*sin(x)+sin(2x)+21*sin(300x)
x0(1)=2;
x0(2)=1;
x0(300)=21;
% evaluating f at all sample points xx
% f(1), f(4), f(6).....f(340) f(356) f(400)
xx=[1 4 6 12 54 69 75 80 89 132 133 152 178 230 300 340 356 400];
% C is measurement matrix
for i=1:n
C(:,i)=sin(i*xx)';
end
b = C*x0;

% b is the result of evaluating f at all sample points xx
% f(1), f(4), f(6).....f(340) f(356) f(400)
% C is the measurement matrix
% let us solve for x and see if it is close to x0
% solve l1-minimization using SeDuMi

cvx_begin
    variable x(n);
    minimize(norm(x,1));
    subject to
        C*x == b;
cvx_end

figure(1)
plot(abs(x-x0),'o')
title('Error between solution solution found by L1 and actual solution')
figure(2)
plot(x,'*')
hold
plot(x0,'o')
title('Solution found for x using L1 minimization')