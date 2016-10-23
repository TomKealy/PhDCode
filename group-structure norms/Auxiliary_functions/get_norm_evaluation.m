function [ norm_u ] = get_norm_evaluation( u, G, D )

%GET_NORM_EVALUATION
%
% get_norm_evaluation( u, G, D )
%
% INPUTS :
%
%          u  : p-dimensional vector 
%          G  : group (logical) matrice (size p x ng, with ng the number of groups)
%               G(i,j) = true if variable i is in group j, false otherwise
%          D  : weight matrice (same size as G)
%               D(i,j) = weight of the variable i in the group j if i is in
%               j, 0 otherwise
%
% OUTPUS :
%
%          norm_u : evaluate the norm omega of vector u
%
%   Copyright (c) 2009 Rodolphe Jenatton. All rights reserved.

ng = size(G,2);

norm_u = 0;

for g=1:ng,
    
    norm_u = norm_u + norm( u( G(:,g) ) .* D( G(:,g), g ), 2);
    
end



end
