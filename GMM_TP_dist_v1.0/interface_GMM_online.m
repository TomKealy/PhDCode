function Xpre = interface_GMM_online(para)
load('data_global');
[para.row, para.col, para.T] = size(C);
para.M = size(Y,3);
para.patchSize = [8 8 8]; 

patchSize = para.patchSize;
n = prod(patchSize);
delta1 = patchSize/2;
delta2 = patchSize/2;
para.R = 1e-6;
para.C = 20;



A = Phi2patches_fast(C, para.patchSize(1), para.patchSize(2), para.T, delta1, delta2);
load('GMM_model_C20_P8X8X8');

model.Mu = Mu;
model.Sig = Sig;
model.pai = pai;

Xt = []; para.kappa = 2;
for m = 1:para.M
    y = video2patches_fast(Y(:,:,m), para.patchSize(1), para.patchSize(2), 1, delta1, delta2);
    X = GMM_3D_JBY(y, A, para, model);
    save('temp','X');
    load('temp');
    Xt = [Xt X];
    Xpre(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(X, para.row, para.col, para.patchSize(1), para.patchSize(1), para.T, delta1, delta2);
    if m < para.M
        [pai, Mu, Sig] = GMM_online_update_kappa(Xt,para,size(X,2),m);
        model.Mu = Mu;
        model.Sig = Sig;
        model.pai = pai;
    end
end







