function [X,W,Fobj,time]=dicolearningAL(signal,TailleDic,lambdasparsecode,lambdadico,mu,mudico,options)


% Dictionary learning using Augmented Lagrangian method.
%
%
%
% options.nbitermax     : number of alternate optimization iterations
% options.threshold     : stopping criterion for the Dico Learning of wrt
%                         relative variation of dictionary Frobenius norm
% options.eta          : factor of reduction of the AL threshold
%                           over the DL iterations
% options.warmstart     : do warm-start over the DL iterations
%
%
% options.fista.nbitermax  : number of iteration in the AL algorithm
% options.fista.threshold  : stopping criterion wrt relative variation of
%                         objective value
% options.fista.thresholdmin : minimal stopping criterion
%

% A. Rakotomamonjy 25/12/10
dim=size(signal,1);
if ~isfield(options,'Xinit');
    Xi=rand(dim,TailleDic);
    Xi=Xi./(ones(dim,1)*sqrt(sum(Xi.^2)));
    X=Xi;
else
    X=options.Xinit;
end;
%mudico=1;

% check if there is any TV regularizer
if strcmp(options.matapprox.leftregularizer,'admmvect')
    for i=1:length(options.matapprox.leftregparam.operatorsadmm);
        if strcmp(options.matapprox.leftregparam.operatorsadmm(i).type,'tv')
            tvregularizer=1;
           lambdatv= options.matapprox.leftregparam.operatorsadmm(i).parameters.bound;
           mudico=lambdatv*10;
           break
        else
            tvregularizer=0;
            lambdatv=0;
            %mudico=0;
        end;
    end;
else
    tvregularizer=0;
    lambdatv=0;
end

tic
mini=inf;
optionsauglag.threshold=options.matapprox.threshold;
optionsauglag.rightregularizer=options.matapprox.rightregularizer;
optionsauglag.leftregularizer=options.matapprox.leftregularizer;
optionsauglag.nbitermax=options.matapprox.nbitermax;
for i=1:options.nbitermax
    optionsauglag.threshold=optionsauglag.threshold*options.eta;
    optionsauglag.threshold=max([optionsauglag.threshold options.matapprox.thresholdmin]);
    optionsauglag.direction='right';
    optionsauglag.regularizerparam= options.matapprox.rightregparam;

    if i>1 & options.warmstart
        optionsauglag.init=W;%+ randn(size(W))*0.0 ;
    end;
    [W,F]=MatrixApproxAL(X,signal,lambdasparsecode,mu,optionsauglag);

    if sum(sum(abs(W)))==0;
        W=randn(TailleDic,size(signal,2));
    end;
    optionsauglag.direction='left';
    optionsauglag.regularizerparam= options.matapprox.leftregparam;
    if i>1 & options.warmstart
        optionsauglag.init=X;% +randn(size(X))*0.0;
    end;
    Xold=X;
    X=MatrixApproxAL(W,signal,lambdadico,mudico,optionsauglag);
    time(i)=toc;     
    
    Fobj(i)=0.5*norm(signal-X*W,'fro').^2+lambdasparsecode*(sum(sum(abs(W))));
        if lambdatv>0
        Fobj(i)=Fobj(i)+lambdatv*sum(sum(abs(diff(X))));
    end;
    Fobj(i)
    if norm(Xold-X,'fro')/norm(Xold,'fro')<options.threshold
        break
    end;  
%  if i>=2 & Fobj(i)>Fobj(i-1)
  %      X=Xold;
    %        optionsauglag.direction='right';
      %      [W,F]=MatrixApproxAL(X,signal,lambdasparsecode,mu,optionsauglag);
       % break
   % end;
    
%     if size(X,2)~=TailleDic
%         Xi=randn(dim,TailleDic-size(X,2));
%         Xi=Xi./(ones(dim,1)*sqrt(sum(Xi.^2)));
%         X=[X Xi];
%     end;

end;

