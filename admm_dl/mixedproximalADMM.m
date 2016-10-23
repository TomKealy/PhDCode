function [u]=mixedproximalADMM(v,operators,options)


% operators(i).type
% operators(i).parameters

% options.nbitermax

nbconstraint=length(operators);
mu=1;
z=zeros(size(v,1),nbconstraint);
y=zeros(size(v,1),nbconstraint);
u=zeros(size(v));
for i=1:options.nbitermax

    % step 1 of ADMM
    uold=u;
    u=1/(nbconstraint*mu+1)*(v+mu*sum(y+z,2));
    if norm(u-uold)< options.seuil
        break
    end;
    for j=1:nbconstraint
        switch operators(j).type
            case 'unit-norm'
                if norm(u-z(:,j))>1
                    y(:,j)=(u-z(:,j))/norm(u-z(:,j)); % unit-norm
                else
                    y(:,j)=(u-z(:,j));
                end;
                %z(:,j)=z(:,j)+y(:,j)-u;
            case 'tv'
                boundtv=operators(j).parameters.bound;
                nbitertv=operators(j).parameters.nbitertv;
                thresholdtv=operators(j).parameters.threshold;
                %   optionstv=operators(j).parameters.optionstv;
                %                 optionstv.u=y(:,j);
                %                 aux=u-z(:,j);
                %                 optionstv.u=y(:,j);
                %y(:,j)=perform_tv_projection_fb(aux,boundtv,optionstv);

                %z(:,j)=z(:,j)+y(:,j)-u;
                aux=u-z(:,j);

                if i==1 || ~operators(j).parameters.warmstart
                 [y(:,j),divpinit,p1,p2] = tvdenoise(aux,mu/boundtv,nbitertv,thresholdtv);
                 else
                  [y(:,j),divpinit,p1,p2] = tvdenoise(aux,mu/boundtv,nbitertv,thresholdtv,divpinit,p1,p2);
                 end;



         case 'tv2'
                sizeim=12;
                boundtv=operators(j).parameters.bound;
                nbitertv=operators(j).parameters.nbitertv;
                thresholdtv=operators(j).parameters.threshold;
                %   optionstv=operators(j).parameters.optionstv;
                %                 optionstv.u=y(:,j);
                %                 aux=u-z(:,j);
                %                 optionstv.u=y(:,j);
                %y(:,j)=perform_tv_projection_fb(aux,boundtv,optionstv);

                %z(:,j)=z(:,j)+y(:,j)-u;
                aux=u-z(:,j);

                if i==1 || ~operators(j).parameters.warmstart
                [auxr,divpinit,p1,p2] = tvdenoise(reshape(aux,sizeim,sizeim),mu/boundtv,nbitertv,thresholdtv);
                y(:,j)=reshape(auxr,sizeim*sizeim,1);
                else
                 [auxr,divpinit,p1,p2] = tvdenoise(reshape(aux,sizeim,sizeim),mu/boundtv,nbitertv,thresholdtv,divpinit,p1,p2);
                                 y(:,j)=reshape(auxr,sizeim*sizeim,1);
                end;
                
            case 'sparse_unit-norm'
                thresh=operators(j).parameters.threshold/mu;
                aux=u-z(:,j);
                y(:,j)= sign(aux).*max(0,abs(aux)-thresh);
                if norm(y(:,j))>1
                    y(:,j)= y(:,j)/norm(y(:,j));
                end;
                
            case 'sparse'
                          thresh=operators(j).parameters.threshold/mu;
                aux=u-z(:,j);
                y(:,j)= sign(aux).*max(0,abs(aux)-thresh);
            case 'positive'
                y(:,j)=((u-z(:,j))>0).*(u-z(:,j)); % prox positive
        end;
        z(:,j)=z(:,j)+y(:,j)-u;
    end;
end;
