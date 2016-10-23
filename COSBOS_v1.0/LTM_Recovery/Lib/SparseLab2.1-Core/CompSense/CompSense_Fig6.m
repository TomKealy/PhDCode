%================================================
% This function produces the results of Figure 6 in the paper, 
% testing the results for varying number of measurements,
% before and after the optimization of the projections
%================================================

function []=CompSense_Fig6()

% Setting Parameters
K=120; % number of atoms in the dictionary
n=80; % dimension of the signals
q=100000; % number of examples
qT=300; % number of error to accumulate
s=4; % cardinality of the sparse representation
m=16:2:40; % number of projections

% Creation of the dictionary
D=randn(n,K);
D=D*diag(1./sqrt(diag(D'*D)));

% Creation of test signals having sparse representations
A=zeros(K,q); A=sparse(A);
h=waitbar(0,'Building the sparse representations ...');
for j=1:1:q,
    waitbar(j/q);
    pos=randperm(K);
    pos=pos(1:s);
    A(pos,j)=randn(s,1);
end;
close(h);
X=D*A;

RESULT=zeros(length(m),4);
for jjj=1:1:length(m), % using different number of projections

    % Generating the compressed measurements
    P=randn(m(jjj),n);
    for j=1:1:m(jjj),
        P(j,:)=P(j,:)/norm(P(j,:));
    end;
    Y=P*X; % the measurements

    Dhat=P*D;
    out=zeros(q,1);
    h=waitbar(0,'Experimenting ...');
    for j=1:1:q,
        waitbar(max(j/q,sum(out(1:j))/qT));
        out(j)=BP_Test(Dhat,Y(:,j),A(:,j),2);
        if out(j)>=0.9999,
            out(j)=BP_Test(Dhat,Y(:,j),A(:,j),1);
        else,
            out(j)=0;
        end;
        if sum(out(1:j))>=qT,
            break;
        end;
    end;
    RESULT(jjj,1)=mean(out(1:j));
    close(h);

    %  decoding using OMP
    Ahat=A*0;
    Dhat=P*D;
    out=zeros(q,1);
    for j=1:1:q,
        Ahat(:,j)=OMP(Dhat,Y(:,j),s);
        out(j)=mean((Ahat(:,j)-A(:,j)).^2)/mean(A(:,j).^2)>1e-3;
        if sum(out(1:j))>=qT,
            break;
        end;
    end;
    RESULT(jjj,2)=mean(out(1:j));

    % optimizing the projection
    iter=1000; dd1=0.8; dd2=0.95;
    Pnew=DesignProjection(D,m(jjj),iter,dd1,dd2,P);

    Dhat=Pnew*D; Y=Pnew*X;
    out=zeros(q,1);
    h=waitbar(0,'Experimenting ...');
    for j=1:1:q,
        waitbar(max(j/q,sum(out(1:j))/qT));
        out(j)=BP_Test(Dhat,Y(:,j),A(:,j),2);
        if out(j)>=0.9999,
            out(j)=BP_Test(Dhat,Y(:,j),A(:,j),1);
        else,
            out(j)=0;
        end;
        if sum(out(1:j))>=qT,
            break;
        end;
    end;
    RESULT(jjj,3)=mean(out(1:j));
    close(h);
    
    % the greedy new results
    Ahat=A*0;
    Dhat=Pnew*D;
    Y=Pnew*X;
    out=zeros(q,1);
    for j=1:1:q,
        Ahat(:,j)=OMP(Dhat,Y(:,j),s);
        out(j)=mean((Ahat(:,j)-A(:,j)).^2)/mean(A(:,j).^2)>1e-3;
        if sum(out(1:j))>=qT,
            break;
        end;
    end;
    RESULT(jjj,4)=mean(out(1:j));

    figure(1); clf; 
    Y1=RESULT(:,1); Y2=RESULT(:,3);
    semilogy(m,Y1,'b'); hold on; semilogy(m,Y2,'b--');
    Y1=RESULT(:,2); Y2=RESULT(:,4);
    semilogy(m,Y1,'r'); hold on; semilogy(m,Y2,'r--');
    drawnow;

end;

Y1=RESULT(:,1); Y2=RESULT(:,3);
semilogy(m,Y1,'b'); hold on; semilogy(m,Y2,'b:');
Y1=RESULT(:,2); Y2=RESULT(:,4);
semilogy(m,Y1,'b'); hold on; semilogy(m,Y2,'b--');
xlabel('m');
ylabel('Relative # of errors');
grid on;
%gtext('OMP with random P');
%gtext('OMP with optimized P');
%gtext('BP with ');
%gtext('random P');
%gtext('BP with optimized P');

%=======================================================

function [OUT]=BP_Test(R,b,x,method)

[k,n]=size(R);
niter=150;
mu0=1e-10;
mumax=1e-5;

if method==1,
    xt=linprog(ones(2*n,1),[],[],[R, -R],b,zeros(2*n,1),ones(2*n,1)*100);
    xt=xt(1:n)-xt(n+1:end);
    OUT=mean((xt-x).^2)/mean(x.^2)>1e-3;
else,
    J=find(abs(x)>1e-10);
    JJ=find(abs(x)<=1e-10);
    nj=length(J);
    njj=length(JJ);
    RJ=R(:,J);
    RJJ=R(:,JJ);
    RR=[RJ RJJ];
    s=sign(x(J));
    z=zeros(size(JJ));
    mu=mu0*ones(njj,1);
    for iter=1:niter
        u=[RJ mulmd(RJJ,mu)]'\[s;z];
        RJJtu=RJJ' * u;
        mu= min( max(1.05*mu.*abs(RJJtu),mu0), mumax);
        if max(abs(RJJtu))<1.000001, break; end
    end
    OUT=max(abs(RJJtu));
end;

%=======================================================

function a=OMP(D,x,L)

% Orthonormal Matching Pursuit with L non-zeros 

[n,K]=size(D);
a=[];
residual=x;
indx=zeros(L,1);
for j=1:1:L,
    proj=D'*residual;
    pos=find(abs(proj)==max(abs(proj)));
    pos=pos(1);
    indx(j)=pos;
    a=pinv(D(:,indx(1:j)))*x;
    residual=x-D(:,indx(1:j))*a;
end;
temp=zeros(K,1); 
temp(indx)=a;
a=sparse(temp);

return;

%========================================================

function [P]=DesignProjection(D,n,Iter,dd1,dd2,Init)
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
%
% Outputs:  P - The projection matrix
%
% Example: N=100; L=200;
%               D=randn(N,L);
%               D=D*diag(1./sqrt(diag(D'*D))); 
%               n=20; Iter=100; dd1=0.7; dd2=0.95; 
%               P=randn(n,N);
%               P=DesignProjection(D,n,Iter,dd1,dd2,P); 
%========================================

[N,L]=size(D);
disp(['The best achievable \mu is ',num2str(sqrt((L-n)/(n*(L-1))))]);
% pause; 

if nargin==5,
    P=randn(n,N); % initialization
    % P=P*diag(1./sqrt(diag(P'*P))); % normalize columns
else, 
    P=Init;
end;
G=D'*P'*P*D; % compute the Gram matrix of the projected dictionary
G=diag(1./sqrt(diag(G)))*G*diag(1./sqrt(diag(G))); % nromalize columns
gg=sort(abs(G(:))); 
pos=find(abs(G(:))>gg(round(dd1*(L^2-L))) & abs(G(:)-1)>1e-6);
RefVal=mean(abs(G(pos)));

for k=1:1:Iter,
    
    % shrink the high inner products
    gg=sort(abs(G(:))); 
    pos=find(abs(G(:))>gg(round(dd1*(L^2-L))) & abs(G(:)-1)>1e-6);
    G(pos)=G(pos)*dd2;
    
    % expand the near zero products
    % pos=find(abs(G(:))<gg(round((1-dd1)*(L^2-L))));
    % G(pos)=G(pos)/dd2;

    % reduce the rank back to n
    method=2;
    if method==1, 
        [U,S,V]=svds(G,n); 
        S=S(1:n,1:n);
        U=U(:,1:n);
    elseif method==2,
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
    pos=find(abs(G(:))>gg(round(dd1*(L^2-L))) & abs(G(:)-1)>1e-6);
    if RefVal>mean(abs(G(pos))), 
        RefVal=mean(abs(G(pos)));
    else,
        return;
    end;
    fprintf(1,'%6i  %12.8f  %12.8f \n',[k,mean(abs(G(pos))),max(abs(G(pos)))]);

end;

return;






