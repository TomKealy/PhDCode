function [X,W,Fobj,itertot]=dicolearningmajorization(signal,TailleDic,lambdasparsecode,options)

% Dictionary learning using majorization method.
%
% This is an implementation of the Yaghoobi's algorithm.
% It actually uses a simple iterative thresholding algorithms and
% a proximal operator.
%
%
% options.nbitermax     : number of alternate optimization iterations
% options.threshold     : stopping criterion for the Dico Learning of wrt
%                         relative variation of dictionary Frobenius norm
%
% options.it.nbitermax  : number od iteration in the IT algorithm
% options.it.threshold  : stopping criterion wrt relative variation of
%                         objective value
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
if ~isfield(options,'Winit');
    W=zeros(TailleDic,size(signal,2));
else
    W=options.Winit;
end;
itertot=0;
Foldobj=inf;
for i=1:options.nbitermax
    
    
    
    
    Wold=W+1;
    Xold=X+1;
    pas=1/norm(X'*X);
    k=1;
    Fold=0.5*norm(signal-X*W,'fro').^2;
    F=0;
    while k<options.it.nbitermax & (abs(Fold-F)/abs(Fold))>options.it.threshold;
        abs(Fold-F);
        Wold=W;
        W= W + pas*X'*(signal-X*W);
        W=sign(W).*(abs(W)-lambdasparsecode*pas).*((abs(W)-lambdasparsecode*pas)>0);
        Fold=F;
        F=0.5*norm(signal-X*W,'fro').^2 + lambdasparsecode*sum(sum(abs(W)));;
        k=k+1;
        if sum(sum(abs(W)))==0;
            W=randn(TailleDic,size(signal,2));
            break
        end;
    end;
    itertot=itertot+k-1;
    pas=1/norm(W*W');
    k=1;
    Fold=norm(signal-X*W,'fro').^2*0.5;
    F=0;
    Xold=X;
    while k<options.it.nbitermax & (abs(Fold-F)/abs(Fold))>options.it.threshold;
        X=X+ pas *(signal-X*W)*W';
        denom=sqrt(sum(X.^2));
        ind=find(denom>0);
        X(:,ind)=X(:,ind)./(ones(size(X,1),1)*denom(ind));
        Fold=F;
        F=0.5*norm(signal-X*W,'fro').^2 ;
        k=k+1;
    end;
    itertot=itertot+k-1;
      Fobj(i)=0.5*norm(signal-X*W,'fro').^2+lambdasparsecode*(sum(sum(abs(W))));
   %   Fobj(i)
    %varnorm(i)=norm(Xold-X,'fro')/norm(Xold,'fro');
    
%     if norm(Xold-X,'fro')/norm(Xold,'fro')<options.threshold
%         break
%     end;
        if (abs(Foldobj-Fobj(i))/abs(Foldobj))<options.threshold;
        break
    end;  
    Foldobj=Fobj(i);

    
end;
% if isempty(Fobj)
%     Fobj(i)=0.5*norm(signal-X*W,'fro').^2+lambdasparsecode*(sum(sum(abs(W))));
% end;
