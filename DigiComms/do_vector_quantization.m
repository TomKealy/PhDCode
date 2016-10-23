%do_vector_quantization.m
clear, clf
var = 0.1;  sdev = sqrt(var); 
Ms=[100 100 100 100 100 100 100 100];
randn('state',0)
class_centers =[1 1;-1 1;-1 -1;1 -1];
N = 2; % Dimension of the sample vectors
K = 4; % Number of clusters
X = [];
for k=1:K  % To create a sample matrix
   X = [X; kron(ones(Ms(k),1),class_centers(k,:))+sdev*randn(Ms(k),N)];
end
Tol = 1e-12; % Tolerance on the quantization error
[codebook,class,quant_errs] = ??????_????????????(?,?,???); 
