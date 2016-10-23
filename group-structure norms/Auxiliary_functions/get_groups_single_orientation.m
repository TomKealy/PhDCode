function [ G, D, IDX ] = get_groups_single_orientation( Height, Width, params )

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



%==========================================================================
%==========================================================================

%--------------------------------------------------------------------------
% Process params
%--------------------------------------------------------------------------
if isfield( params, 'slope' ),
    slope = abs(params.slope);
else
    slope = 0;
end
%--------------------------------------------------------------------------
if isfield( params, 'exp_factor' ),
    exp_factor = params.exp_factor;
else
    exp_factor = 0.5;
end
%--------------------------------------------------------------------------
if isfield( params, 'cardG_weights' ),
    cardG_weights = params.cardG_weights;
else
    cardG_weights = false;
end
%--------------------------------------------------------------------------
if isfield( params, 'intersection' ),
    intersection = params.intersection;
else
    intersection = false;
end
%==========================================================================
%==========================================================================

ncells = Height * Width;

[ rows, cols ] = ndgrid( Height:-1:1, 1:Width );

cols = cols(:);
rows = rows(:);

if isinf( slope ) || ( slope == 0 ),    
    EQ1 = cols;
    EQ2 = rows;
else
    EQ1 = rows - slope * cols;
    EQ2 = rows + slope * cols;
end

%==========================================================================
%==========================================================================
% Groups from the top-left corner to the bottom-right corner
% i.e., rows - slope * cols >= constant (we have to remove the smallest constant not to have the full group, i.e., POSITION(1)  )

POSITION = unique( EQ1 ); % Note that unique sorts the result
L        = length(POSITION);

G1 = false( ncells, L );
D1 = zeros( ncells, L );

for j=1:L,

    index       = ( EQ1 >= POSITION(j) );
    G1(index,j) = true;
    
    if cardG_weights, 
        
        D1(index,j) = ncells / sum(index)^2 ;%D1(index,j-1) = ncells / sum(index)^2 ;
        
    else

        for i=j:L,

            D1( ( EQ1 == POSITION(i) ) , j) = exp_factor^(i-j);

        end

    end

end

%--------------------------------------------------------------------------
% Groups from the the bottom-right corner to the top-left one
% i.e., rows - slope * cols <= constant (we have to remove the largest constant not to have the full group, i.e., POSITION(end) )

G2 = false( ncells, L );
D2 = zeros( ncells, L );

for j=1:L,

    index       = ( EQ1 <= POSITION(j) );
    G2(index,j) = true;
    
    if cardG_weights, 
        
        D2(index,j) = ncells / sum(index)^2 ;
        
    else
        
        for i=1:j,

            D2( ( EQ1 == POSITION(i) ) , j) = exp_factor^(j-i);

        end

    end
    
end

%==========================================================================
% Groups from the the bottom-right corner to the top-left one
% i.e., rows + slope*cols <= constant (we have to remove the largest constant not to have the full group, i.e., POSITION(end) )

POSITION = unique( EQ2 ); % Note that unique sorts the result
L        = length(POSITION);

G3 = false( ncells, L );
D3 = zeros( ncells, L );

for j=1:L,

    index       = ( EQ2 <= POSITION(j) );
    G3(index,j) = true;
    
    if cardG_weights, 
        
        D3(index,j) = ncells / sum(index)^2 ;
        
    else
        
        for i=1:j,

            D3( ( EQ2 == POSITION(i) ) ,j) = exp_factor^(j-i);

        end

    end
    
end

%--------------------------------------------------------------------------
% Groups from the the top-right corner to the bottom-left corner
% i.e., rows + slope*cols >= constant (we have to remove the smallest constant not to have the full group, i.e., POSITION(1)  )

G4 = false( ncells, L );
D4 = zeros( ncells, L );

for j=1:L,

    index       = ( EQ2 >= POSITION(j) );
    G4(index,j) = true;
    
    if cardG_weights, 
        
        D4(index,j) = ncells / sum(index)^2 ;
        
    else
        
        for i=j:L,

            D4( ( EQ2 == POSITION(i) ) , j) = exp_factor^(i-j);

        end

    end

end

%==========================================================================
%==========================================================================
if cardG_weights,
    D1 = D1 / max( D1(:) );
    D2 = D2 / max( D2(:) );  
    D3 = D3 / max( D3(:) );    
    D4 = D4 / max( D4(:) ); 
end;
%==========================================================================
%==========================================================================
% If we do not perform islasso, we have to remove the full group for each orientation

if ~intersection,    

  G1 = G1(:,2:end);
  D1 = D1(:,2:end);
  
  G2 = G2(:,1:end-1);
  D2 = D2(:,1:end-1);
     
  G3 = G3(:,1:end-1);
  D3 = D3(:,1:end-1);
  
  G4 = G4(:,2:end);
  D4 = D4(:,2:end);
 
end
%==========================================================================
%==========================================================================
if ( Height == 1 ), 
    G3 = []; G4 = []; D3 = []; D4 = [];  
end;
if ( Width  == 1 ), 
    G1 = []; G2 = []; D1 = []; D2 = []; 
end;
%==========================================================================
%==========================================================================

[trash, sort_idx1] = sort( sum(G1,1), 'descend' );

[trash, sort_idx2] = sort( sum(G2,1), 'descend' );

[trash, sort_idx3] = sort( sum(G3,1), 'descend' );

[trash, sort_idx4] = sort( sum(G4,1), 'descend' );


G = logical( [G1(:,sort_idx1) G2(:,sort_idx2) G3(:,sort_idx3) G4(:,sort_idx4)] );
D = [D1(:,sort_idx1) D2(:,sort_idx2) D3(:,sort_idx3) D4(:,sort_idx4)];

%==========================================================================
%==========================================================================

IDX = []; m = 0;
if ~isempty(G1),
    IDX = [IDX ; [1,size(G1,2)] ];
    m   = size(G1,2);
end
if ~isempty(G2), 
    IDX = [IDX ; [1,size(G2,2)]+m ]; 
    m   = m+size(G2,2); 
end
if ~isempty(G3), 
    IDX = [IDX ; [1,size(G3,2)]+m ]; 
    m = m+size(G3,2); 
end
if ~isempty(G4), 
    IDX = [IDX ; [1,size(G4,2)]+m ];
end

end

