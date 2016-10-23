function [M_f,x_f,y_f,alpha]=l1_ls_random_data_producer(m,n,L,M_mean,M_var,x_mean,x_div,h_mean,h_var,Dist,zero_per,lambda_vec,L_f,max_iter,seed)
%   M_F(L*m,n)     M_f(m,n,L)
%   x_F(n,1)       x_f(n,L)
%   y_F(L*m,1)     y_f(m,L)
%   h_F(L*m,1)     h_f(m,L) 
init_rand(seed);
temp_V=random('Exponential',1,n,1);
temp_Z=random('Normal',0,1,n,1);
temp_X=x_mean+x_div*sqrt(2*temp_V).*temp_Z;
[~,temp_ind]=sort(abs(temp_X),'ascend');
temp_X(temp_ind(1:floor(zero_per*n)))=0;
x_F=(Dist/norm(temp_X,2))*temp_X;
x_f=x_F*ones(1,L); %#ok<NASGU>

M_f=random('Normal',M_mean,M_var,m,n,L);
M_f_T(1:n,1:m,1:L)=0;
for l=1:L
    M_f_T(1:n,1:m,l)=(sqrt(L_f)/norm(M_f(1:m,1:n,l),2))*M_f(1:m,1:n,l)';
    M_f(1:m,1:n,l)=M_f_T(1:n,1:m,l)';
end;
M_F=(reshape(M_f_T,n,L*m))';

h_F=random('Normal',h_mean,h_var,L*m,1);
h_f=reshape(h_F,m,L); %#ok<NASGU>

y_F=M_F*x_F+(eye(m*L)-M_F*pinv(M_F))*h_F; 
y_f=reshape(y_F,m,L);

alpha=1/norm(M_F,2)^2;

x_CS=x_F;

for k=1:max_iter
    x_CS=wthresh(x_CS-alpha*M_F'*(M_F*x_CS-y_F),'s',alpha*lambda_vec);
end;

if norm(x_CS-x_F,2)/n>0.1
    warning('LASSO might not give the true recovery');
    disp(strcat('norm(x_CS-x_F,2)=',num2str(norm(x_CS-x_F,2)/n)));
    
end;

plot(x_F);
hold all;
plot(x_CS);
legend('x_F','x_{CS}');
% keyboard;
close all;

x_f=x_CS*ones(1,L);

end


