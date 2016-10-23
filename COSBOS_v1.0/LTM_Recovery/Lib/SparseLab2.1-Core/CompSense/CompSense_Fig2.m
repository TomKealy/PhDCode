%================================================
% This function generates Figure 2 in the paper, whcih shows the 
% shrink operation that needs to be performed.
%================================================

function []=CompSense_Fig2()

t=0.5; gamma=0.6;
g=-1:0.001:1;
h=(abs(g)>=t).*(g*gamma)+...
    (abs(g)<t & abs(g)>=gamma*t).*(sign(g)*gamma*t)...
    +(abs(g)<gamma*t).*g;
plot(g,g,':');
grid on;
axis([-1 1 -1 1])
hold on;
plot(g,gamma*g,'b:')
axis image
h=plot(g,h);
set(h,'LineWidth',2); 
xlabel('Input Value');
ylabel('Output Value');