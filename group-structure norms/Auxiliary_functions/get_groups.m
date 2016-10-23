function [ G, D, IDX ] = get_groups( Height, Width, params )

% GENERATE GROUP & WEIGHT MATRICES FOR 2-D GRIDS
%
% The 2-d grid is encoded as 
%             
%       |1       height+1  2*height+1 ... (width-1)*height+1|
%       |2                            ...                   |
%       |.                            ...                   |
%       |.                            ...                   |
%       |.                            ...                   |
%       |height  2*height             ...       width*height|
%
%
% FOR AN EXAMPLE OF APPLICATION, SEE THE Demo_for_groups.m FILE.
%
%
% get_groups( Height, Width, params )
%
% INPUTS:
%
% Height      : height of the grid
% Width       : width of the grid
% params      : structure that contains additional paramaters (see the fields below)
%
%       params.slope: orientation of the groups (by default, 0)
%                     If ( abs(params.slope) == Inf | abs(params.slope) ==  0), we have the rectangular groups
%                     If (abs(params.slope) == 1), we have the +- Pi/4 groups
%
%                     If length(params.slope)>1, then the group/weight matrices G and D are built for the
%                     length(params.slope) orientations
%  
%       params.exp_factor: weight parameter (by default, 0.5; see the beginning of the experiment section in the paper,
%                          such weights are denoted by (W3). Note that if params.exp_factor=1, we have the
%                          uniform weights, denoted by (W1) in the paper) 
%
%       params.cardG_weights: weight parameter (by default, false)
%                     (see the beginning of the experiment section in the paper, such weights are denoted by (W2))
%
%       params.intersection: weight parameter (by default, false)
%                     needs to be put to true for ISlasso
%
%
% OUTPUTS:
% 
% G           : group (logical) matrix (size p x ng, with ng the number of groups)
%               G(i,j) = true if variable i is in group j, false otherwise
% D           : weight matrix (same size as G)
%               D(i,j) = weight of the variable i in the group j if i is in j, 0 otherwise
% IDX         : index matrix that sums up the sorting of G, per direction and size
%               (size nd x 2, with nd the number of directions in G; for
%               instance 4 for the rectangular groups)
%
%               For example:
%                   IDX(d,1):IDX(d,2) correspond to the column indexes of groups in G for the direction d. 
%     
%   Copyright (c) 2009 Rodolphe Jenatton. All rights reserved.

G   = logical([]);
D   = [];
IDX = [];

%--------------------------------------------------------------------------
if isfield( params, 'slope' ),
    
    for slope=params.slope,
        
        params_        = params;
        params_.slope  = slope;
        
        [ G_, D_, IDX_ ] = get_groups_single_orientation( Height, Width, params_ );

       
        G   = [G, G_];%#ok<AGROW>
        D   = [D, D_];%#ok<AGROW>
        
        if isempty(IDX),
            IDX = IDX_;
        else
            IDX = [IDX; max(IDX(:))+IDX_]; %#ok<AGROW>
        end
        
    end
             
else
    
    [ G, D, IDX ] = get_groups_single_orientation( Height, Width, params );
    
end
%--------------------------------------------------------------------------


end
