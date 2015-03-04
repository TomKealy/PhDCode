function Xpre = GMM_3D_JBY(Y, CS, para,model)

Mu = model.Mu;
Sig = model.Sig;
pai = model.pai;

[d,Nt] = size(Y);
Ri = para.R*eye(d);


% GMM Inversion
Xpre = zeros(size(CS{1},2),Nt);

% for i = 1:Nt
parfor i = 1:Nt  % (parallel computing for multicore)
    if mod(i,1000) == 1
        fprintf('i=%d\n',i);
    end
    [pai1,temp,~] = GMM_Inference(Y(:,i), CS{i}, Ri, Mu, Sig, pai);
    Xpre(:,i) = temp;
end


