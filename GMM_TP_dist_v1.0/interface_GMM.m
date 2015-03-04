function Xpre = interface_GMM(para)
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

for m = 1:para.M
    y = video2patches_fast(Y(:,:,m), para.patchSize(1), para.patchSize(2), 1, delta1, delta2);
    X = GMM_3D_JBY(y, A, para, model);
    Xpre(:,:,(m-1)*para.T+1:m*para.T) = patches2video_fast(X, para.row, para.col, para.patchSize(1), para.patchSize(1), para.T, delta1, delta2);   
end







