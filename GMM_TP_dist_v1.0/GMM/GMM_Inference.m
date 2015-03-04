function [paii,X_bar, mu_M] = GMM_Inference(Y,Phi,Ri,mu,Sig,pai)        
[p,c] = size(mu); 
evd = zeros(1,c); 
mu_M = zeros(p,c);
m = length(Y); 
for t = 1:c
    mu1 = mu(:,t);
    Q = Sig(:,:,t)*Phi';
    P = inv(Phi*Q+Ri);
%     P = eye(m)/(Phi*Q+Ri);
    P = (P+P')/2;
    Sig1 = Q*P;
    res = Y - Phi*mu1;
    tmp = -0.5*m*log(2*pi) + sum(log(diag(chol(P)))) - 0.5*sum(res.*(P*res),1);

    mu_M(:,t) = mu1 + Sig1*res;
    evd(:,t) = tmp + log(pai(t)+eps);
end
[paii,evm] = Design_Pai(evd); 
X_bar = mu_M*paii';


function [paii,evm] = Design_Pai(evd)
c = size(evd,2); 
mmax  = max(evd,[],2); 
paii = exp(evd+repmat(-mmax,[1,c]))+eps; 
evm = sum(paii,2); 
paii = paii./repmat(evm,[1,c]); 
evm = log(evm) + mmax;


