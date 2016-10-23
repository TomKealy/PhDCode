function  [W,F,V,D]=MatrixApproxAL(X,signal,lambda,mu,options)
%[W,F,V,D]=MatrixApproxAL(X,signal,lambda,mu,options)
% Augmented Lagragian for Matrix approximation
%
% direction = right
% min_{W}  1/2 || S- XW||_F^2 + \lambda \Omega(W)
%
% direction = left
% min_{W}  1/2 || S- WX||_F^2 + \lambda \Omega(W)
%
%

if isfield(options,'initW');
    W=options.initW;
    V=options.initV;
    D=options.initD;
else

    switch options.direction
        case 'right'
            K=size(signal,2);
            TailleDic=size(X,2);
            W=zeros(TailleDic,K);
            options.regularizer=options.rightregularizer;
            V=zeros(size(W));
            D=zeros(size(W));
        case 'left'
            K=size(X,1);
            W=ones(size(signal,1),K);
            denom=sqrt(sum(W.^2));
            ind=find(denom>0);
            W(:,ind)=W(:,ind)./(ones(size(W,1),1)*denom(ind));
            %pas=1/norm(X*X');
            V=zeros(size(W));
            D=zeros(size(V));
            options.regularizer=options.leftregularizer;
    end;
end;

% adjust regularizer parameters  in some proximal operator
lambdasparse=0;
lambdatv=0;
if strcmp(options.regularizer,'admmvect')
    for i=1:length(options.regularizerparam.operatorsadmm);
        if strcmp(options.regularizerparam.operatorsadmm(i).type,'sparse') || ...
                strcmp(options.regularizerparam.operatorsadmm(i).type,'sparse_unit-norm')
            options.regularizerparam.operatorsadmm(i).parameters.threshold= ...
                options.regularizerparam.operatorsadmm(i).parameters.threshold/mu;
            lambdasparse=options.regularizerparam.operatorsadmm(i).parameters.threshold;

        end;
        if strcmp(options.regularizerparam.operatorsadmm(i).type,'tv')
            options.regularizerparam.operatorsadmm(i).parameters.bound= ...
                options.regularizerparam.operatorsadmm(i).parameters.bound*mu;
            lambdatv=options.regularizerparam.operatorsadmm(i).parameters.bound;
        end;
        if strcmp(options.regularizerparam.operatorsadmm(i).type,'tv2')
            options.regularizerparam.operatorsadmm(i).parameters.bound= ...
                options.regularizerparam.operatorsadmm(i).parameters.bound*mu;
            lambdatv=options.regularizerparam.operatorsadmm(i).parameters.bound;
        end;
    end;
end

%V=W;
Wold=W+1;

az=1;


switch options.direction
    case 'right'
        invXtX=inv(X'*X+mu*eye(size(X,2)));
        R=invXtX*(X'*signal);
        Fold=0.5*norm(signal-X*W,'fro').^2;
    case 'left'
        invXXt=inv(X*X'+mu*eye(size(X,1)));
        R=(signal*X')*invXXt;
        %  XXt=X*X';
        %  SXt=signal*X';
        %  pas=1/norm(XXt);
        Fold=0.5*norm(signal-W*X,'fro').^2;
end;





F=[];
j=0;
while j <= 2 || (abs(F(j-1)-F(j))/abs(F(j-1)))>options.threshold
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
    %  First alternate step
    %----------------------------------------------------------------------



    W=V+D;
    switch options.direction
        case 'right'
            W=R+mu*invXtX*W;
        case 'left'
            W=R+mu*W*invXXt;
    end;

    thresh=lambda/mu;
    V=W-D;
    %end;
    if j > options.nbitermax
        break
    end;
    %----------------------------------------------------------------------
    % Second alternate  proximal step
    %----------------------------------------------------------------------


    switch options.regularizer
        case '\ell_1-\ell_2'
            K=size(V,2);
            normeV=sqrt(sum(V.^2,2));      % shrink
            bool=normeV>(thresh);
            aux=(1-thresh./normeV).*bool;
            V=V.*(aux*ones(1,K));           % shrink
        case '\ell_1-\ell_1'
            V=sign(V).*(abs(V)-thresh).*((abs(V)-thresh)>0);
        case 'frobenius2'
            V=V/(1+2*lambda);
        case 'proj-unitnorm2'       % projection on the l2 unit-ball
            denom=sqrt(sum(V.^2));
            ind=find(denom>1);
            V(:,ind)=V(:,ind)./(ones(size(V,1),1)*denom(ind));

        case 'unit-norm_tv';
            type ='unit-norm_tv';
            parameters = options.regularizerparam;
            VV=mixedproximal(V,type,parameters);
            V=VV;
        case 'admmvect'
            for kk=1:size(V,2);
                VV(:,kk)=mixedproximalADMM(V(:,kk),options.regularizerparam.operatorsadmm,options.regularizerparam.optionsadmm);
            end;
            V=VV;
        case 'none'
    end;

    D=D+V-W;



    j=j+1;


end;

