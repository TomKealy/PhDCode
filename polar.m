x = [101101];
y = zeros(1,length(x));
n=5
for i=1:length(x)
     if(x(i)==1)
         y(i*n:(i+1)*n)=1
     else
         y(i*n:(i+1)*n)=-1
     end
end
stairs(y)