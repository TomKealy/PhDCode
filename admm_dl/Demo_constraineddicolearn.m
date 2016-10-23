%
%   Demo code for constrained dictionary learning
%
%



clear all
close all
%
addpath('./ksvd/');
addpath('./ksvd/');
addpath('./ksvd/ompbox/');
addpath('./ksvd/tools/');
dim=100;
TailleDic=10;
T=3; % nb de fonctions de bases
K=200; % nb de signaux
nbitermax=1;
rsnr=30;
verbose=1;


rand('state',0);
randn('state',0);



for iter=1:nbitermax
    
    %--------------------------------------------------------------------------
    % creating dictionary and signals
    %--------------------------------------------------------------------------
    
    
    Xt=randn(dim,TailleDic);
    
    i=1;
    while i <= TailleDic
        %for i=1:TailleDic
        Xt(:,i)=zeros(dim,1);
        deb=floor(rand*dim/2);
        long=floor(rand*dim);
        Xt(deb+1:min(deb+long,dim),i)=cos((deb+1:min(deb+long,dim))/(rand*20) + rand*pi);
        if norm(Xt(:,i))>0
            i=i+1;
        end;
    end;
    %Xt=abs(Xt);
    
    Xt=Xt./(ones(dim,1)*sqrt(sum(Xt.^2)));
    
    
    
    poidsM=zeros(TailleDic,K);
    for i=1:K
        ind=randperm(size(Xt,2));
        indice=ind(1:T);
        poids=rand(T,1)*0.8+0.2;
        signalt(:,i)=Xt(:,indice)*poids;
        poidsM(indice,i)=poids;
    end;
    
    stdnoise=std(signalt)/rsnr;
    signal=signalt+randn(size(signalt)).*(ones(dim,1)*stdnoise);
    
    %----------------------------------------------------------------------
    %       Options
    %----------------------------------------------------------------------
    lambdasparsecode=0.05;lambdadico=0;
    options.nbitermax=50;
    options.threshold=1e-4;    
    options.warmstart=0;
    options.eta=1;
    
    % options maj
    options.it.threshold=1e-8;
    options.it.nbitermax=200000;
    
    % options fista/AL
    options.matapprox.nbitermax=1000;
    options.matapprox.threshold=1e-8;  
    options.matapprox.thresholdmin=1e-8
    options.matapprox.rightregularizer='\ell_1-\ell_1';
    options.matapprox.rightregparam=0;
    
%-------------------------------------------------------------------------    
%   Projection on unit-norm regularizer
%-------------------------------------------------------------------------
%     options.matapprox.leftregularizer='proj-unitnorm2';
%     options.matapprox.leftregparam=0;
    
%--------------------------------------------------------------------------
%   TV regularization  + sparse unit-norm  
%--------------------------------------------------------------------------
    options.matapprox.leftregularizer='admmvect';
    operators(1).type='tv';% a diminuer si moins de bruit
    operators(1).parameters.nbitertv=300;
    operators(1).parameters.threshold=1e-3;
    operators(1).parameters.warmstart=1;
    operators(1).parameters.bound=0.01;

    operators(2).type='sparse_unit-norm';
    operators(2).parameters.threshold=0.01;   % a diminuer si moins de bruit
   % operators(3).type='positive';
    

    optionsadmm.nbitermax=500;
    optionsadmm.seuil=1e-3;
    options.matapprox.leftregparam.optionsadmm=optionsadmm;
    options.matapprox.leftregparam.operatorsadmm=operators;
    
    
    %----------------------------------------------------------------------
    %           Majorization method
    %----------------------------------------------------------------------
    %    here a sparse constraint on the weight and  a 
    %    tic

    %     [Xm,Wm,Fobjw,itertot]=dicolearningmajorization(signal,TailleDic,lambdasparsecode,options);
    %
    %     timemaj(iter)=toc;
    %     [distm(iter),ratiom(iter)]=dictdist(Xm,Xt);
    %
    
    
    
    
    %----------------------------------------------------------------------
    %           FISTA method
    %----------------------------------------------------------------------
%     tic
%     [Xf,Wf,Fobj,time]=dicolearningFista(signal,TailleDic,lambdasparsecode,lambdadico,options);
%     timefista(iter)=toc;
%     [distfista(iter),ratiof(iter)]=dictdist(Xf,Xt)
%     
    
    
    %----------------------------------------------------------------------
    %           AL method
    %----------------------------------------------------------------------
    mu=1;
    mudico=10;
    tic
    [Xal,Wal,Fobjal,time]=dicolearningAL(signal,TailleDic,lambdasparsecode,lambdadico,mu,mudico,options);
    timeal(iter)=toc;
    [distal(iter),ratio]=dictdist(Xal,Xt)
end;

plot(Xal);
title('Learned dictionary atoms');

