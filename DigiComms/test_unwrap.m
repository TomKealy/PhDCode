%test_unwrap.m
% What is the function of unwrap()?
n=[0:200]; th=2*pi*n/100; x=exp(j*th);
th1=angle(x);  th2=unwrap(th1);
plot(th,th1,'b:', th,th2,'r-')
