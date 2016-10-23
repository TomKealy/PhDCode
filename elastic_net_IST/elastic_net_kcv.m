function [la_kcv,err_kcv] = elastic_net_kcv(X,Y,eps,range,k,split,err_type)
%%%%%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%%%%%%
%parameter choice through K-fold or LOO cross validation for elastic net
%regularization
%
%   [LA_KCV] = ELASTIC_NET_KCV(X,Y,EPS,RANGE) returns the regularization
%   parameter which minimizes LOO cross-validation of elastic net 
%   regularization with the training set (X, Y) with relative parameter 
%   EPS. If RANGE is a 2x1 vector LA_KCV will be chosen among the values in
%   the geometric series RANGE(1).^RANGE(2), if RANGE is a 3x1 
%   vector LA_KCV will be chosen among the RANGE(3) values in the geometric 
%   series from RANGE(1) to RANGE(2), if RANGE is a 1x3 vector LA_KCV will 
%   be chosen among the values in the equispaced series 
%   RANGE(1):RANGE(2):RANGE(3),otherwise LA_KCV will be chosen among the 
%   values in RANGE. 
% 
%   [LA_KCV, ERR_KCV] = ELASTIC_NET_KCV(X,Y,EPS,RANGE) also returns the 
%   cross-validation error
% 
%   [...] = ELASTIC_NET_KCV(X,Y,EPS,RANGE,K) performs sequential K-fold 
%   cross-validation. If K=0 or K=length(Y) it performs LOO cross-validation
% 
%   [...] = ELASTIC_NET_KCV(X,Y,EPS,RANGE,K,SPLIT) performs 
%   sequential(SPLIT=0) K-fold cross-validation, or random (SPLIT=1) K-fold 
%   cross-validation.
%
%   [...] = ELASTIC_NET_KCV(X,Y,EPS,RANGE,K,SPLIT,ERR_TYPE) select the type
%   of error to be considered for KCV. If ERR_TYPE= 'CLASS' classification 
%   error is considered (default), if ERR_TYPE='REGR' the squared error is 
%   considered. 
%
if nargin<4, error('too few input'), end
if nargin<5, k = 0; end; 
if nargin<6, split = 0; end; 
if nargin<7, err_type = 'class'; end
if nargin>7, error('too many input'), end

% evaluates range of values for the parameter
if isequal(size(range),[2,1]); % geometric series I
    lambda =1/sqrt(length(Y))*range(1).^(1:range(2)); 
elseif isequal(size(range),[3,1]); % geometric series II
   lambda = range(1)*((range(2)/range(1))^(1/range(3))).^(1:range(3)); 
elseif isequal(size(range),[1,3]); % equispaced series
    lambda = range(1):range(2):range(3); 
else lambda = range;
end

sets = splitting(Y,k,split); %splits the training set in k subsets

err_kcv = zeros(length(lambda),1); % initialization

for i = 1:length(sets);
    fprintf('\n split number %d',i);
    ind = setdiff(1:length(Y),sets{i}); % indexes of the i-th training set
    Xtr = X(ind,:); % input matrix of the i-th training set
    Ytr = Y(ind); % output vector of the i-th training set
    Xts = X(sets{i},:); % input matrix of the i-th test set
    Yts = Y(sets{i}); % output vector of the i-th test set
    % for each value of the parameter, evaluates test error
    for l = 1:length(lambda);
        beta = elastic_net_lrn(Xtr,Ytr,eps,lambda(l));
        err_kcv(l) = err_kcv(l)+elastic_net_test(Xts,Yts,beta,err_type)/k;
    end
end

% determines the optimal parametr as the largest one minimizing the
% cross-validation error
la_kcv = lambda(find(err_kcv==min(err_kcv),1,'last'));

function sets = splitting(Y,K,type)
%SPLITTING Splitting in balanced subsets
%   SETS = SPLIT(Y,K) return a cell array of K subsets of 1:n where
%   n=length(Y). The elements 1:n are split so that in each subset the
%   ratio between indexes corresponding to positive elements of array Y and
%   indexes corresponding to negative elements of Y is the about same as in
%   1:n. The subsets are obtained  by sequentially distributing the
%   elements of 1:n.
%
%   SETS = SPLIT(Y,K,TYPE) return a cell array of K subsets of 1:n where
%   n=length(Y) according to splitting type TYPE which can be either
%   0 (sequential) or 1 (random). SPLIT(Y,K,0) = SPLIT(Y,K)
%
if nargin==2, type = 0; end; 

n = length(Y);
if or(K==0,K==n);
    sets = cell(1,n);
    for i = 1:n, sets{i} = i; end
else
    c1 = find(Y>=0);
    c2 = find(Y<0);
    l1 = length(c1);
    l2 = length(c2);
    if type==0;
        perm1=1:l1;
        perm2=1:l2;
    elseif type==1;
        perm1 = randperm(l1);
        perm2 = randperm(l2);
    end;
    sets = cell(1,K);
    i = 1;
    while i<=l1;
        for v = 1:K;
            if i<=l1;
                sets{v} = [sets{v}; c1(perm1(i))];
                i = i+1;
            end;
        end;
    end;
    i = 1;
    while i<=l2;
        for v = 1:K;
            if i<=l2;
                sets{v} = [sets{v}; c2(perm2(i))];
                i = i+1;
            end;
        end;
    end;
end    