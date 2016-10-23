function [X,W,Fobj,time]=dicolearningFista(signal,TailleDic,lambdasparsecode,lambdadico,options)

% Dictionary learning using FISTA method.
%
%
%
% options.nbitermax     : number of alternate optimization iterations
% options.threshold     : stopping criterion for the Dico Learning of wrt
%                         relative variation of dictionary Frobenius norm
% options.eta          : factor of reduction of the FISTA threshold
%                           over the DL iterations
% options.warmstart     : do warm-start over the DL iterations
%
%
% options.matapprox.nbitermax  : number of iteration in the FISTA algorithm
% options.matapprox.threshold  : stopping criterion wrt relative variation of
%                         objective value
% options.matapprox.thresholdmin : minimal stopping criterion
%

% A. Rakotomamonjy 25/11/10


dim=size(signal,1);
if ~isfield(options,'Xinit');
    Xi=rand(dim,TailleDic);
    Xi=Xi./(ones(dim,1)*sqrt(sum(Xi.^2)));
    X=Xi;
else
    X=options.Xinit;
end;

% check if there is any TV regularizer
if strcmp(options.matapprox.leftregularizer,'admmvect')
    for i=1:length(options.matapprox.leftregparam.operatorsadmm);
        if strcmp(options.matapprox.leftregparam.operatorsadmm(i).type,'tv')
            tvregularizer=1;
           lambdatv= options.matapprox.leftregparam.operatorsadmm(i).parameters.bound;
           break
        else
            tvregularizer=0;
            lambdatv=0;
        end;
    end;
else
    tvregularizer=0;
    lambdatv=0;
end


tic
optionsfista.threshold=options.matapprox.threshold;
optionsfista.rightregularizer=options.matapprox.rightregularizer;
optionsfista.leftregularizer=options.matapprox.leftregularizer;
optionsfista.nbitermax=options.matapprox.nbitermax;
if isfield(options.matapprox,'accelerated');
    optionsfista.accelerated=options.matapprox.accelerated;
end;
for i=1:options.nbitermax
    optionsfista.threshold=optionsfista.threshold*options.eta;
    optionsfista.threshold=max([optionsfista.threshold options.matapprox.thresholdmin]);
    optionsfista.direction='right';
    optionsfista.regularizerparam= options.matapprox.rightregparam;
    %
    if i>1 & options.warmstart
        optionsfista.init=W;%+ randn(size(W))*0.0 ;
    end;
    optionsfista.direction='right';
    [W,F]=MatrixApproxFista(X,signal,lambdasparsecode,optionsfista);
    
 
   

    optionsfista.direction='left';
    optionsfista.regularizerparam= options.matapprox.leftregparam;
    if i>1 & options.warmstart
        optionsfista.init=X;% +randn(size(X))*0.0;
    end;
    Xold=X;
%     indact=find(sum(abs(W),2)>0);
%     W=W(indact,:);
    %Xold=X(:,indact);
    X=MatrixApproxFista(W,signal,lambdadico,optionsfista);
    time(i)=toc; 
    Fobj(i)=0.5*norm(signal-X*W,'fro').^2+lambdasparsecode*(sum(sum(abs(W))));
    if lambdatv>0
        Fobj(i)=Fobj(i)+lambdatv*sum(sum(abs(diff(X))));
    end;
    Fobj(i)
    if norm(Xold-X,'fro')/norm(Xold,'fro')<options.threshold
        break
    end;
 %   if i>=2 & Fobj(i)>Fobj(i-1)
  %      X=Xold;
    %        optionsfista.direction='right';
     %       [W,F]=MatrixApproxFista(X,signal,lambdasparsecode,optionsfista);
     %   break
  %  end;
    
%     if length(indact)~=TailleDic
%         Xi=randn(dim,TailleDic-length(indact));
%         Xi=Xi./(ones(dim,1)*sqrt(sum(Xi.^2)));
%         X=[X Xi];
%     end;



end;
