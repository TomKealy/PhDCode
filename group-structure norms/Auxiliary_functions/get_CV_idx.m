function [ idx_train, idx_test ] = get_CV_idx( N, KFolds )

%GET_CV_IDX
%   Copyright (c) 2009 Rodolphe Jenatton. All rights reserved.

rd_idx = randperm(N);
L      = N/KFolds;


idx       = cell(KFolds,1);
idx_train = cell(KFolds,1);
idx_test  = cell(KFolds,1);



for k=1:KFolds,
    idx{k} = rd_idx( floor((k-1)*L)+1:floor(k*L) );
end

for k=1:KFolds,
    idx_test{k} = idx{k};
    for j=1:(k-1), 
        idx_train{k} =  [idx_train{k} idx{j}];
    end
    for j=(k+1):KFolds, 
        idx_train{k} =  [idx_train{k} idx{j}];
    end
end