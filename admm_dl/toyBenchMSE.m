function toyBenchMSE(comment,learning,options);

% % This is the script that has been used for comparing the MSE
% % efficiency of KSVD, MAJ ISTA, and FISTA.
% % it has to be launcher with launcherMSE.m



% clear all
% close all
% optionsdefault;





dim=learning.dim;        % problem dimension
TailleDic=learning.TailleDic;   % Dictionary size
T=learning.T;            % nb de fonctions de bases
K=learning.K;           % nb de signaux
nbitermax=learning.nbitermax;
rsnr=learning.rsnr;
verbose=learning.verbose;
type=learning.type;
nbtest=learning.nbtest;
nbval=learning.nbval;
lambdasparsecodevec=learning.lambdasparsecodevec; % regularization parameters
lambdadico=learning.lambdadico;


filename=['MSE-' comment 'Dictype-' type '-dim' int2str(dim) '-TailleDic' ...
    int2str(TailleDic) '-nbsig' int2str(K) '-rsnr' int2str(learning.rsnr) ...
    '-mu' int2str(options.ALmu) '-mudico' int2str(options.ALmudico)]
remtime(filename,0,'Initialisation')

if isfield(options,'comment')
    if isfield(options.comment,'tv')
        filename=[filename '-' options.comment.tv];
    end;
    if isfield(options.comment,'sparse')
        filename=[filename '-' options.comment.sparse];
    end;
    filename=[filename '.mat'];
end;
nbitermaxaux=options.nbitermax;
for iter=1:nbitermax
    iter
    rand('state',iter);
    randn('state',iter);
    %----------------------------------------------------------------------
    % creating dictionary and signals
    %----------------------------------------------------------------------

    switch type
        case 'gaussian'
            Xt=randn(dim,TailleDic);

        case 'smooth'
            Xt=randn(dim,TailleDic);
            i=1;
            while i <= TailleDic
                %for i=1:TailleDic
                Xt(:,i)=zeros(dim,1);
                deb=floor(rand*dim*0.75);
                long=floor(rand*dim/3)+dim/5;
                Xt(deb+1:min(deb+long,dim),i)=cos((deb+1:min(deb+long,dim))/(rand*20) + rand*pi);
                if norm(Xt(:,i))>0
                    i=i+1;
                end;
            end;
            %Xt=abs(Xt);
            Xt=Xt./(ones(dim,1)*sqrt(sum(Xt.^2)));
        case 'constant'
            i=1;
            while i <= TailleDic
                length=floor(rand*dim/3)+dim/6;
                deb=floor(rand*(dim-length));
                Xt(:,i)=zeros(dim,1);
                Xt(deb+1:deb+length,i)=rand+1;
                if norm(Xt(:,i))>0
                    i=i+1;
                end;
            end;
    end;
    Xt=Xt./(ones(dim,1)*sqrt(sum(Xt.^2)));

    nbdata=K+nbtest+nbval;
    poidsM=zeros(TailleDic,nbdata);
    for i=1:nbdata
        ind=randperm(size(Xt,2));
        indice=ind(1:T);
        poids=rand(T,1)*0.8+0.2;
        signalt(:,i)=Xt(:,indice)*poids;
        poidsM(indice,i)=poids;
    end;

    stdnoise=std(signalt)/rsnr;
    signal=signalt+randn(size(signalt)).*(ones(dim,1)*stdnoise);
    indperm=randperm(nbdata);
    indperm=1:nbdata;
    signalapp=signal(:,indperm(1:K));
    signalval=signal(:,indperm(K+1:K+nbval));
    signaltest=signal(:,indperm(K+nbval+1:nbdata));
    signaltesttrue=signalt(:,indperm(K+nbval+1:nbdata));
    % dictionary initialization
    Xi=rand(dim,TailleDic);
    Xi=Xi./(ones(dim,1)*sqrt(sum(Xi.^2)));
    options.Xinit=Xi;


    %----------------------------------------------------------------------
    %           Majorization method
    %----------------------------------------------------------------------
    nbitermaxaux=options.nbitermax;
    options.nbitermax=options.nbitermax*10;
    mini=inf;
    for iterlambda=1:size(lambdasparsecodevec,2);
        lambdasparsecode=lambdasparsecodevec(iterlambda);
        tic
        [Xm,Wm,Fobjw{iter,iterlambda},timew{iter,iterlambda}]=dicolearningmajorization(signalapp,TailleDic,lambdasparsecode,options);
        timemaj(iter,iterlambda)=toc;
        [distm(iter,iterlambda),ratiom(iter,iterlambda)]=dictdist(Xm,Xt);
        options.matapprox.direction='right';
        [W,F]=MatrixApproxFista(Xm,signalval,lambdasparsecode,options.matapprox);
        msemval(iterlambda)=0.5*norm(Xm*W-signalval,'fro').^2;
        if msemval(iterlambda)<mini
            Xmopt=Xm;
            lambdaopt=lambdasparsecode;
            mini= msemval(iterlambda);
        end;
    end;
    distmopt(iter)=dictdistance(Xmopt,Xt);
    [W,F]=MatrixApproxFista(Xmopt,signaltest,lambdaopt,options.matapprox);
    msem(iter)=0.5/nbtest*norm(Xmopt*W-signaltesttrue,'fro').^2;
    lambdaoptmaj(iter)=lambdaopt;
    Xmmat{iter}=Xmopt;
    
    
    % % % % %
    % % % %
    % % % %     %----------------------------------------------------------------------
    % % % %     %           Majorization method
    % % % %     %----------------------------------------------------------------------
    %     nbitermaxaux=options.nbitermax;
    %     options.nbitermax=options.nbitermax*10;
    %     mini=inf;
    %     for iterlambda=1:size(lambdasparsecodevec,2);
    %         lambdasparsecode=lambdasparsecodevec(iterlambda);
    %         tic
    %         [Xm,Wm,Fobjw2{iter,iterlambda},timew2{iter,iterlambda}]=dicolearningmajorization(signalapp,TailleDic,lambdasparsecode,options);
    %         timemaj2(iter,iterlambda)=toc;
    %         [distm2(iter,iterlambda),ratiom(iter,iterlambda)]=dictdist(Xm,Xt);
    %         options.matapprox.direction='right';
    %         [W,F]=MatrixApproxFista(Xm,signalval,lambdasparsecode,options.matapprox);
    %         msemval2(iterlambda)=0.5*norm(Xm*W-signalval,'fro').^2;
    %         if msemval2(iterlambda)<mini
    %             Xmopt2=Xm;
    %             lambdaopt2=lambdasparsecode;
    %             mini= msemval(iterlambda);
    %         end;
    %     end;
    %     distmopt2(iter)=dictdistance(Xmopt,Xt);
    %     [W,F]=MatrixApproxFista(Xmopt2,signaltest,lambdaopt2,options.matapprox);
    %     msem2(iter)=0.5/nbtest*norm(Xmopt2*W-signaltesttrue,'fro').^2;
    %     lambdaoptmaj2(iter)=lambdaopt;
    %
    %  [W]=(Xmopt2'*Xmopt2)\(Xmopt2'*signaltest);
    % msem2(iter)=0.5/nbtest*norm(Xmopt2*W-signaltesttrue,'fro').^2;


    % % % %     %----------------------------------------------------------------------
    % % % %     %           KSVD method
    % % % %     %----------------------------------------------------------------------

    params.data = signalapp;
    params.Tdata = T;
    params.dictsize = TailleDic;
    params.iternum = 1000;
    params.memusage = 'high';
    params.exact=0;


    mini=inf;
    Tvec=[2:10];
    for iterT=1:size(Tvec,2);
        tic
        params.Tdata = Tvec(iterT);
        [Dksvd,g,err] = ksvd(params,'');
        timeksvd(iter,iterT)=toc;
        [distksvd(iter,iterT),ratiom(iter,iterT)]=dictdist(Dksvd,Xt);
        for iterlambda=1:size(lambdasparsecodevec,2);
            lambdasparsecode=lambdasparsecodevec(iterlambda);
            options.matapprox.direction='right';
            [W,F]=MatrixApproxFista(Dksvd,signalval,lambdasparsecode,options.matapprox);
            mseksvdval(iterT,iterlambda)=0.5*norm(Dksvd*W-signalval,'fro').^2;
            if mseksvdval(iterT,iterlambda)<mini
                Dksvdopt=Dksvd;
                Topt2=Tvec(iterT);
                lambdaopt=lambdasparsecode;
                mini= mseksvdval(iterT);
            end;
        end;
    end;
    distdksvdopt(iter)=dictdistance(Dksvd,Xt);
    [W,F]=MatrixApproxFista(Dksvdopt,signaltest,lambdaopt,options.matapprox);
    mseksvd(iter)=0.5/nbtest*norm(Dksvdopt*W-signaltesttrue,'fro').^2;
    lambdaoptksvd(iter)=lambdaopt;
    Topt(iter)=Topt2;
Xksvdmat{iter}=Dksvdopt;


    % % % %     %----------------------------------------------------------------------
    % % % %     %           FISTA method
    % % % %     %----------------------------------------------------------------------
    %       mini=inf;
    %     options.nbitermax=nbitermaxaux;
    %     for iterlambda=1:size(lambdasparsecodevec,2)
    %         lambdasparsecode=lambdasparsecodevec(iterlambda);
    %         tic
    %         [Xf,Wf,Fobj{iter,iterlambda},timef{iter,iterlambda}]=dicolearningFista(signalapp,TailleDic,lambdasparsecode,lambdadico,options);
    %         timefista(iter,iterlambda)=toc;
    %         objfista(iter,iterlambda)=Fobj{iter,iterlambda}(end);
    %         [distfista(iter,iterlambda),ratiof(iter,iterlambda)]=dictdist(Xf,Xt);
    %         options.matapprox.direction='right';
    %         [W,F]=MatrixApproxFista(Xf,signalval,lambdasparsecode,options.matapprox);
    %         msefval(iterlambda)=0.5*norm(Xf*W-signalval,'fro').^2;
    %         if msefval(iterlambda)<mini
    %             Xfopt=Xf;
    %             lambdaopt=lambdasparsecode;
    %             mini= msefval(iterlambda);
    %         end;
    %     end;
    %     distfistaopt(iter)=dictdistance(Xfopt,Xt);
    %     [W,F]=MatrixApproxFista(Xfopt,signaltest,lambdaopt,options.matapprox);
    %     msefista(iter)=0.5/nbtest*norm(Xfopt*W-signaltesttrue,'fro').^2;
    %    lambdaoptfista(iter)=lambdaopt;
    %      Xfmat{iter}=Xfopt;
    %
    %     %----------------------------------------------------------------------
    %     %           ISTA method
    %     %----------------------------------------------------------------------
    %     mini=inf;
    %     options.nbitermax=nbitermaxaux;
    %     options.matapprox.accelerated=0;
    %     for iterlambda=1:size(lambdasparsecodevec,2)
    %         lambdasparsecode=lambdasparsecodevec(iterlambda);
    %             options.matapprox.accelerated=0;
    %         tic
    %
    %         [Xf,Wf,Fobji{iter,iterlambda},timei{iter,iterlambda}]=dicolearningFista(signalapp,TailleDic,lambdasparsecode,lambdadico,options);
    %         timeista(iter,iterlambda)=toc;
    %        objista(iter,iterlambda)=Fobji{iter,iterlambda}(end);
    %         [distista(iter,iterlambda),ratiof(iter,iterlambda)]=dictdist(Xf,Xt);
    %         options.matapprox.direction='right';
    %             options.matapprox.accelerated=1;
    %         [W,F]=MatrixApproxFista(Xf,signalval,lambdasparsecode,options.matapprox);
    %         msefival(iterlambda)=0.5*norm(Xf*W-signalval,'fro').^2;
    %         if msefival(iterlambda)<mini
    %             Xfopt=Xf;
    %             lambdaopt=lambdasparsecode;
    %             mini= msefival(iterlambda);
    %         end;
    %     end;
    %     dististaopt(iter)=dictdistance(Xfopt,Xt);
    %     [Wi,F]=MatrixApproxFista(Xfopt,signaltest,lambdaopt,options.matapprox);
    %     mseista(iter)=0.5/nbtest*norm(Xfopt*Wi-signaltesttrue,'fro').^2;
    %    lambdaoptista(iter)=lambdaopt;
    %     Ximat{iter}=Xfopt;
    % %
    %----------------------------------------------------------------------
    %           AL method
    %----------------------------------------------------------------------
    mini=inf;
    options.nbitermax=nbitermaxaux;
    for iterlambda=1:size(lambdasparsecodevec,2)
        lambdasparsecode=lambdasparsecodevec(iterlambda);
        tic
        [Xal,Wal,Fobjal{iter,iterlambda},timeal{iter,iterlambda}]=dicolearningAL(signalapp,TailleDic,lambdasparsecode,lambdadico,options.ALmu,options.ALmudico,options);
        timead(iter,iterlambda)=toc;
        objal(iter,iterlambda)=Fobjal{iter,iterlambda}(end);
        [distal(iter,iterlambda),ratio]=dictdist(Xal,Xt);
        [W,F]=MatrixApproxFista(Xal,signalval,lambdasparsecode,options.matapprox);
        msealval(iterlambda)=0.5*norm(Xal*W-signalval,'fro').^2;
        if msealval(iterlambda)<mini
            Xalopt=Xal;
            lambdaopt=lambdasparsecode;
            mini= msealval(iterlambda);
        end;
    end;
    distalopt(iter)=dictdistance(Xalopt,Xt);
    [W,F]=MatrixApproxFista(Xalopt,signaltest,lambdaopt,options.matapprox);
    mseal(iter)=0.5/nbtest*norm(Xalopt*W-signaltesttrue,'fro').^2;
    lambdaoptal(iter)=lambdaopt;
    Xalmat{iter}=Xalopt;
    %
    

    %
     Xtrue{iter}=Xt;
       signalappmat{iter}=signalapp; 
       save(filename)

end;
clear signal signalapp signalval signaltest signalt
save(filename)
%save('truedicosmooth','Xtrue');
%dim

