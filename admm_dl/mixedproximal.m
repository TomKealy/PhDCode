function Y=mixedproximal(X,type,parameters)

% positive sparse unit norm
%
%       type.one =  positive;           x positive
%       type.two = sparse_unit-norm;    \lambda \ell_1 st ||x||_2 = 1
%
%
% unit-norm bounded total variation
%
%       type ='unit-norm_tv'           parameters.bound = bound on TV
%
%
% parameters.nbiterdougrach     : nb of iteration max for douglas-radford
% parameters.threshold : threshold for the DR proximal algorithm

nbitermax=parameters.nbiterdougrach;
seuil=parameters.threshold;

%--------------------------------------------------------------------------
%
%             positive and sparse unit-norm
%
%--------------------------------------------------------------------------

switch type
    case 'positive_sparse_unit-norm'
        %if (strcmp(type.one,'positive') & strcmp(type.two,'sparse_unit-norm'))
        K=size(X,2);

        lambda=parameters.lambda;


        for i=1:K
            v=X(:,i);
            z=randn(size(v));
            z1old=z;
            z1=z1old+1;
            j=0;
            while j<nbitermax & norm(z1-z1old)>seuil;
                thresh=lambda/2;
                aux= (z+v)/2;
                z1old=z1;
                z1=sign(aux).*max(0,abs(aux)-thresh);
                z1=z1/norm(z1);  % sparse_unit-norm prox
                z2=z-z1+ (2*z1-z).*((2*z1-z)>0); % positive projection
                z=z2;
                j=j+1;
            end;

            Y(:,i)=z1;
        end;


    case 'unit-norm_tv'
        %--------------------------------------------------------------------------
        %
        %               unit-norm bounded total variation
        %
        %--------------------------------------------------------------------------


        K=size(X,2);
        tautv =parameters.bound;
        optionstv.niter = parameters.nbitertv;
        optionstv.u=randn(size(X,1),1);

        for i=1:K
            %i%
            v=X(:,i);
            z=v;%randn(size(v));
            z1old=z;
            z1=z1old+1;
            j=0;
            while j<nbitermax &             norm(z1-z1old)/norm(z1)>seuil;
                thresh=tautv;
                aux= (z+v)/2;
                z1old=z1;
                % tv projection
                %[z1] = perform_tv_projection_fb(aux,thresh, optionstv);
                z1=tvdenoise(aux,1/thresh,100);
                % unit norm projection
                z2=z-z1+ (2*z1-z)/norm(2*z1-z);
                z=z2;
                j=j+1;
                norm(z1-z1old);
            end;
            if norm(z1) ==0
                Y(:,i)=z1;
            else
                Y(:,i)=z1/norm(z1);
            end;
        end;
    otherwise
end;

