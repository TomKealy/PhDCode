%================================================
% This Script produces Figure 3 in the paper, showing the way 
% the proposed iterative algorithm for optimizining P converge 
% for various values of \gamma (used here as dd2).
%================================================

function []=CompSense_Fig3()

% Setting Parameters
K=400; % number of atoms in the dictionary
n=200; % dimension of the signals
m=30;

% Creation of the dictionary
D=randn(n,K);
D=D*diag(1./sqrt(diag(D'*D)));

% Creaiton of the initial projection matrix
P=randn(m,n);

% optimizing the projection

iter=1000; 
dd1=0.2; 
RES=zeros(iter+1,5); 
count=1;
for dd2=0.55:0.1:0.95,
    [Pnew,RES(:,count)]=DesignProjection_OUT(D,m,iter,dd1,dd2,P,2);
    count=count+1;
end;

plot(0:1:iter,RES,'b'); 
%gtext('\gamma=0.55'); 
%gtext('\gamma=0.65'); 
%gtext('\gamma=0.75');
%gtext('\gamma=0.85'); 
%gtext('\gamma=0.95');
xlabel('Iteration'); 
ylabel('Value of \mu_t');

%================================================

function [P,OUTPUT]=DesignProjection_OUT(D,n,Iter,dd1,dd2,Init,methodT)

%========================================
% In this program we iteratively build a projection matrix that when 
% applied to the dictionary D, it leads to better coherence.
% 
% Inputs:    D - The dictionary
%               n - Number of projections (rows in P)
%               Iter - Number of iterations
%               dd1 - Relative number of the G entries to shrink
%               dd2 - Shrink factor to use (1/dd2 is used to expand)
%               Init - An initial projection matrix
%               methodT - threshold method to use
%
% Outputs:  P - The projection matrix
%                 OUTPUT - the \mu_t per each iteration
%
% Example: N=100; L=200;
%               D=randn(N,L);
%               D=D*diag(1./sqrt(diag(D'*D))); 
%               n=20; Iter=100; dd1=0.7; dd2=0.95; 
%               P=randn(n,N);
%               P=DesignProjection(D,n,Iter,dd1,dd2,P); 
%========================================

% Mode of operation
methodSVD=2; % SVD is done by the power-method and not directly
% methodT=1; % working with the percentage option
[N,L]=size(D);
% disp(['The best achievable \mu is ',num2str(sqrt((L-n)/(n*(L-1))))]);

if nargin==5,
    P=randn(n,N); % initialization
    % P=P*diag(1./sqrt(diag(P'*P))); % normalize columns
else, 
    P=Init;
end;
G=D'*P'*P*D; % compute the Gram matrix of the projected dictionary
G=diag(1./sqrt(diag(G)))*G*diag(1./sqrt(diag(G))); % nromalize columns
gg=sort(abs(G(:))); 

if methodT==1,
  Threshold=gg(round(dd1*(L^2-L)));
else, 
  Threshold=dd1;
end;
pos=find(abs(G(:))>Threshold & abs(G(:)-1)>1e-6);
RefVal=mean(abs(G(pos)));
OUTPUT=zeros(Iter+1,1);
OUTPUT(1)=RefVal;

for k=1:1:Iter,
    
    % shrink the high inner products
    gg=sort(abs(G(:))); 
    if methodT==1,
      Threshold=gg(round(dd1*(L^2-L)));
    else,
      Threshold=dd1;
    end;
    pos=find(abs(G(:))>Threshold & abs(G(:)-1)>1e-6);
    G(pos)=G(pos)*dd2;
    pos=find(abs(G(:))<=Threshold & abs(G(:))>Threshold*dd2);
    G(pos)=sign(G(pos))*dd2*Threshold;
    
    % reduce the rank back to n
    if methodSVD==1, 
        [U,S,V]=svds(G,n); 
        S=S(1:n,1:n);
        U=U(:,1:n);
    elseif methodSVD==2,
        U=randn(L,n); 
        for jjj=1:1:10, U=G*U; U=orth(U); end;
        S=diag(diag(U'*G*U)); 
    end;
    G=U*S*U';
    PD=S.^0.5*U';
    P=PD*pinv(D); 
    
    % Normalize the columns
    G=D'*P'*P*D; 
    G=diag(1./sqrt(diag(G)))*G*diag(1./sqrt(diag(G)));

    % Show status
    gg=sort(abs(G(:))); 
    if methodT==1,
      Threshold=gg(round(dd1*(L^2-L)));
    else,
      Threshold=dd1;
    end;    
    pos=find(abs(G(:))>Threshold & abs(G(:)-1)>1e-6);
    % if RefVal>mean(abs(G(pos))), 
        RefVal=mean(abs(G(pos)));
        OUTPUT(k+1)=RefVal;
    % else,
        % return;
    % end;
    fprintf(1,'%6i  %12.8f  \n',[k,mean(abs(G(pos)))]); % ,max(abs(G(pos)))]);

end;

return;

