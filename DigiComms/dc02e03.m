%dc02e03.m 
% Plot Fig. 2.3 to check the validity of the CLT
clear, clf
Ns=10000; NB=100; % the numbers of samples and bins
mx=1/2; % the mean of standard uniform noise: Eq.(2.1.29) (a=0, b=1)
sgmx2=1/12; % the variance of standard uniform noise: Eq.(2.1.30)
sgmx=sqrt(sgmx2); 
Ls=[2 10];
for i=1:length(Ls)
   L=Ls(i); % # of iterations 
   rand('twister',5489); %return rand() to its default initial state
   y = sum(rand(L,Ns))/L; % Eq.(2.1.49)
   hist(y,NB), pause
end
hold on
[ns,cs]=hist(y,NB); 
dy=cs(2)-cs(1); % the bin width
my=mx; sgmy2=sgmx2/L; %the average & variance of the sample averages
sgmy=sqrt(sgmy2);
y=my+[-500:500]*(my/500); % the range on the y-axis
fy=exp(-(y-my).^2/2/sgmy2)/sqrt(2*pi)/sgmy;
plot(y,Ns*dy*fy,'r')
