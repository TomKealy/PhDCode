function [V,d] = eig_r(Q,NumOfNonZero)
% calcuate the largest NumOfNonZero eigenvalues and eigenvectors
[V1,D1] = eig(Q);
d1 = diag(D1);
[val,ind] = sort(d1,'descend');
V = V1(:,ind(1:NumOfNonZero));
d = d1(ind(1:NumOfNonZero));


