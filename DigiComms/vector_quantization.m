function [codebook,class,quant_errs]=vector_quantization(X,K,Tol)
%Finds a KxN codebook representing the samples in an MxN data matrix X,
% of which each row represents a sample, by continuing the iteration 
% until quantization error difference between successive iterations<Tol
% Output arguments
%  codebook: A KxN matrix with each row being a code vector
%  class   : An Mx1 vector of code number tos which each sample belongs
%  quant_erros: Quantization error history
% Example:
%   M=200; N=2; X = randn(M,N); K = 3; Tol = 1e-12;
%   [codebook,class,quantization_error] = vector_quantization(X,K,Tol)
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only

[M,N] = size(X); tmp=randperm(M);
codebook = X(tmp(1:K),:); % Initial the codebook arbitrartily
MaxIter = 1000;  quant_errs = [10; 1]; 
for iter=1:MaxIter
   % Classify each element to the nearest codevector
   codebook_new = zeros(K,N); class = zeros(M,1); 
   membership = zeros(M,K); 
   quant_err = 0;
   for m=1:M
      for k=1:K, dist(k) = sum((X(m,:)-codebook(k,:)).^2); end
      [min_dist,k_min] = min(dist);
      km = k_min(ceil(length(k_min)*rand)); % one from multiple minima
      class(m)=km;  % an Mx1 vector of code numbers for each sample 
      membership(m,km)=1; % an MxK matrix showing the membership       codebook_new(km,:) = codebook_new(km,:) + X(m,:);
      quant_err = quant_err + min_dist;
   end
   for k=1:K
      no_of_elements = sum(membership(:,k));
      if no_of_elements>1
        codebook_new(k,:) = codebook_new(k,:)/no_of_elements;
       elseif no_of_elements==0 
        % if the class represented by a codevector has no element
        codebook_new(k,:) = sum(codebook_new(1:k,:))/k; 
      end 
   end
   codebook = codebook_new; quant_err = quant_err/M;
   quant_errs = [quant_errs; quant_err];
   if abs(quant_errs(end)-quant_errs(end-1))<Tol, break; end
end
