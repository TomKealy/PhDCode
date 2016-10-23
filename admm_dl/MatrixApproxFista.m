function  [W,F]=L2prox(X,signal,lambda,options)

%
%
% min_{W}  1/2 || S- XW||_F^2 + \lambda \Omega(W)
% options.direction='right'
% options.rightregularizer=
%
% or
%
% min_{W}  1/2 || S- WX||_F^2 + \lambda \Omega(W)
% options.direction='left'
% options.leftregularizer=

if isfield(options,'init');
    W=options.init;
    switch options.direction
        case 'right'
            K=size(signal,2);
            TailleDic=size(X,2);
            options.regularizer=options.rightregularizer;
        case 'left'
            K=size(X,1);
            options.regularizer=options.leftregularizer;
    end;
else

    switch options.direction
        case 'right'
            K=size(signal,2);
            TailleDic=size(X,2);
            W=zeros(TailleDic,K);
            options.regularizer=options.rightregularizer;
            %pas=1/norm(X'*X);
        case 'left'
            K=size(X,1);
            W=zeros(size(signal,1),K);
            % denom=sqrt(sum(W.^2));
            %ind=find(denom>0);
            %W(:,ind)=W(:,ind)./(ones(size(W,1),1)*denom(ind));
            %pas=1/norm(X*X');
            options.regularizer=options.leftregularizer;
    end;
end;

if ~isfield(options,'accelerated');
    options.accelerated=1;
end;

V=W;
Wold=W+1;
az=1;


switch options.direction
    case 'right'
        XtX=X'*X;
        XtS=X'*signal;
        pas=1/norm(XtX);
        Fold=0.5*norm(signal-X*W,'fro').^2;
    case 'left'
        XXt=X*X';
        SXt=signal*X';
        pas=1/norm(XXt);

        Fold=0.5*norm(signal-W*X,'fro').^2;
end;


% adjust regularizer parameters  in some proximal operator
lambdasparse=0;
lambdatv=0;
if strcmp(options.regularizer,'admmvect')
    for i=1:length(options.regularizerparam.operatorsadmm);
        if strcmp(options.regularizerparam.operatorsadmm(i).type,'sparse') || ...
                strcmp(options.regularizerparam.operatorsadmm(i).type,'sparse_unit-norm')
            options.regularizerparam.operatorsadmm(i).parameters.threshold= ...
                options.regularizerparam.operatorsadmm(i).parameters.threshold*pas;
            lambdasparse=options.regularizerparam.operatorsadmm(i).parameters.threshold;
        end;
        if strcmp(options.regularizerparam.operatorsadmm(i).type,'tv')
            options.regularizerparam.operatorsadmm(i).parameters.bound= ...
                options.regularizerparam.operatorsadmm(i).parameters.bound*pas;
            lambdatv=options.regularizerparam.operatorsadmm(i).parameters.bound;
        end;
        if strcmp(options.regularizerparam.operatorsadmm(i).type,'tv2')
            options.regularizerparam.operatorsadmm(i).parameters.bound= ...
                options.regularizerparam.operatorsadmm(i).parameters.bound*pas;
            lambdatv=options.regularizerparam.operatorsadmm(i).parameters.bound;
        end;
    end;
end




j=0;

F=[];
while j <= 2 || (abs(F(j-1)-F(j))/abs(F(j-1)))>options.threshold
    %while norm(W-Wold)>options.threshold
    %max(max(abs(XtS-XtX*W)))
    %abs(Fold-F)
    %norm(W-Wold)
    %max(max(abs(W-Wold)))
    switch options.regularizer
        case '\ell_1-\ell_2'
            normeW=sqrt(sum(W.^2,2));
        case '\ell_1-\ell_1'
            normeW=sum(sum(abs(W)));
        case 'frobenius2'
            normeW=norm(W,'fro').^2;
        case 'proj-unitnorm2'
            normeW=0;
        case 'none'
            normeW=0;
        case 'admmvect'
            normeW=0;
            if lambdasparse > 0
                normeW=lambdasparse*sum(sum(abs(W)));
            end;
            if lambdatv>0
                normeW=normeW+lambdatv*sum(sum(abs(diff(W))));

            end;
        otherwise
            normeW=0;
    end;
    %----------------------------------------------------------------------
    % objective value
    %----------------------------------------------------------------------
    switch options.direction
        case 'right'
            fobj=0.5*norm(signal-X*W,'fro').^2;
        case 'left'
            fobj=0.5*norm(signal-W*X,'fro').^2;
    end;
    sumnorm = sum(lambda.*normeW);
    F(j+1)=fobj+sumnorm;

    %----------------------------------------------------------------------
    % gradient and backtracking step
    %----------------------------------------------------------------------

    switch options.direction
        case 'right'
            %grad=-X'*(signal-X*V);
            grad=-XtS+XtX*V;
        case 'left'
            %grad=-(signal-V*X)*X';
            grad=-SXt+V*XXt;
    end;
    %
    %   beware of backtracking computing the trace may be expensive
    %     lin=trace((W-V)'*grad);
    %      the following lines are buggy
    %     while F > fobj+sumnorm+lin + 1/2/pas*norm(W-V,'fro').^2;
    %         pas=pas*eta;
    %     end;
    Wold=W;
    %----------------------------------------------------------------------
    % gradient step
    %----------------------------------------------------------------------
    W=V-pas*grad;   % gradient step
    %----------------------------------------------------------------------
    % proximal step
    %----------------------------------------------------------------------
    switch options.regularizer
        case '\ell_1-\ell_2'
            K=size(W,2);
            normeW=sqrt(sum(W.^2,2));      % shrink
            bool=normeW>(lambda*pas);
            aux=(1-lambda*pas./normeW).*bool;
            W=W.*(aux*ones(1,K));           % shrink
        case '\ell_1-\ell_1'
            W=sign(W).*(abs(W)-lambda*pas).*((abs(W)-lambda*pas)>0);
        case 'frobenius2'
            W=W/(1+2*lambda);
        case 'none'

        case 'proj-unitnorm2'
            denom=sqrt(sum(W.^2));
            ind=find(denom>1);
            W(:,ind)=W(:,ind)./(ones(size(W,1),1)*denom(ind));

        case 'unit-norm_tv';
            type ='unit-norm_tv';
            parameters = options.regularizerparam;
            WW=mixedproximal(W,type,parameters);
            W=WW;
        case 'admmvect'
            for kk=1:size(W,2);
                WW(:,kk)=mixedproximalADMM(W(:,kk),options.regularizerparam.operatorsadmm,options.regularizerparam.optionsadmm);
            end;
            W=WW;
    end;

    %% tseng
    %     a1=2/( (j-1)+3);
    %     delta=W-Wold;
    %     V=W+(1-az)/az*a1*delta;
    %     az=a1;

    %% fista
    %     delta=W-Wold;
    %     a1=(1+sqrt(1+4*az^2))/2;
    %     V=W + (az-1)/a1*delta;
    %     az=a1;
    %     j=j+1;
    %
    if options.accelerated ==1
        delta=W-Wold;
        a1=(1+sqrt(1+4*az^2))/2;
        V=W + (az-1)/a1*delta;
        az=a1;
        j=j+1;
    else

        V=W;
        j=j+1;
    end;


    if j > options.nbitermax
        break
    end;
end;

