function [Phi2,idxCS] = Phi2patches_fast_repeat(Phi, n1, n2, n3, T, delta1, delta2)
% This function convert the 3D cube into 2D patches
% Jianbo Yang
% 06/19/2013

if mod(n1,delta1) ~= 0 || mod(n2,delta2) ~= 0 || mod(size(Phi,1),n1) ~= 0 || mod(size(Phi,2),n2) ~= 0
    error('Set parameters correctly!');
end
A = rand(n1,n2); % in case that all elements in Phi are same
A = repmat(A, size(Phi,1)/n1, size(Phi,2)/n2);
A1 = image2patches_fast(A, n1, n2, delta1, delta2);
[d1, d2] = size(A1);

idxCS  = zeros(d2,1);
a0 = A(1:n1,1:n2); b0 = repmat(a0,2,2);
k = 1;
for j = 0:delta2:n2-1
    for i = 0:delta1:n1-1    
        temp  = b0(1+i:i+n1, 1+j:j+n2);
        PhiAll(:,k) = temp(:)';  
        k = k + 1;
    end
end
k = k - 1;
if rank(PhiAll)<k, error('rank should be greater than k'); end
for i = 1:k
    diff = bsxfun(@minus, A1,  PhiAll(:,i));
    idx = sum(abs(diff)) < 1e-6;
    idxCS(idx) = i;
end




M_patches = zeros(d1,d2,T);
for t = 1:T
    M_patches(:,:,t) = image2patches_fast(Phi(:,:,t), n1, n2, delta1, delta2);
end
n = size(M_patches,2);
Phi2 = cell(k,1);
for i = 1:k
    idx = find(idxCS == i);
    idx = idx(1);
    temp = [];
    for t = 1:T
        temp = [temp sparse(diag(M_patches(:,idx,t)))];
    end
    Phi2{i} = kron(eye(n3/T),temp);
end



return;
