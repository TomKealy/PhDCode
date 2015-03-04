function [pai, Mu, Sig] = GMM_online_update_kappa(Xpre,para,Nt,m);

% Update GMM prior
load([para.dirc])
nd = size(data0,2);

for i = 0:m
    xi(i+1) = (i+1)^para.kappa;
end



xi = xi./sum(xi);
ndt = Nt/xi(end);


nm = round(ndt*xi); 
p = randperm(nd); 
data1 = data0(:,p(1:nd-sum(nm(2:end))));

for i = 1:m
    p = randperm(Nt);
    temp = Xpre(:,(i-1)*Nt+1:(i-1)*Nt+Nt);
    if nm(i+1)<Nt
        data1 = [data1 temp(:,p(1:nm(i+1)))];
    else
        data1 = [data1 temp]; 
    end
end



clear data0 Xpre

options = statset('Display','final','MaxIter',200);
obj = gmdistribution.fit(data1',para.C,'Regularize',10^-3,'Options',options);
clear data1
pai = obj.PComponents;
Mu = obj.mu';
Sig = obj.Sigma;