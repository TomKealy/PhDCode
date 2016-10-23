% There are mostly 3 kinds of group weights:
%
%   1) params.exp_factor = 1, we get the uniform weights, denoted by (W1)
%   in the paper.
%
%   2) params.exp_factor < 1, we get the weights denoted by (W3) in the
%   paper, i.e., the weights that take into account overlaps.
%
%   3) params.cardG_weights = true, we get the weights denoted by (W2) in the 
%   paper, i.e., the weights that depend only on the size of the groups.
%
%
% The group orientations (e.g., rectangular groups) is defined through
% params.slope.
%
%   1) ( abs(params.slope) == Inf | abs(params.slope) ==  0), we have the rectangular groups
%
%   2) ( abs(params.slope) == 1), we have the +- Pi/4 groups
%
%   3) length(params.slope)>1, then the group/weight matrices G and D are built for the length(params.slope) orientations
%                              Example: if params.slope = [0,1], then we have the rectangular AND the +- Pi/4 groups

% We consider a 3 x 4  grid
Height = 3;
Width  = 4;

% By default, rectangular groups with weights params.exp_factor=0.5
% (Denoted by (W3) in the paper)
[G,D, IDX] = get_groups( Height, Width, [] );

% We plot the groups with their corresponding weights for the first orientation
orientation = 1;
for j=IDX(orientation,1):IDX(orientation,2),

    reshape( G(:,j), Height, Width )
    reshape( D(:,j), Height, Width )

end

% We plot the groups with their corresponding weights for the second orientation
orientation = 2;
for j=IDX(orientation,1):IDX(orientation,2),

    reshape( G(:,j), Height, Width )
    reshape( D(:,j), Height, Width )

end

%==========================================================================
%Pi/4 groups with weights denoted by (W3) in the paper.

clear params
params.slope = 1;

[G,D, IDX] = get_groups( Height, Width, params );

% We plot the groups with their corresponding weights for the first orientation
orientation = 1;
for j=IDX(orientation,1):IDX(orientation,2),

    reshape( G(:,j), Height, Width )
    reshape( D(:,j), Height, Width )

end

% We plot the groups with their corresponding weights for the second orientation
orientation = 2;
for j=IDX(orientation,1):IDX(orientation,2),

    reshape( G(:,j), Height, Width )
    reshape( D(:,j), Height, Width )

end
%==========================================================================
%Rectangular groups and Pi/4 groups with weights denoted by (W3) in the paper.

clear params
params.slope = [0,1];

[G,D, IDX] = get_groups( Height, Width, params );

% We plot the groups with their corresponding weights for the first orientation
orientation = 1;
for j=IDX(orientation,1):IDX(orientation,2),

    reshape( G(:,j), Height, Width )
    reshape( D(:,j), Height, Width )

end

% We plot the groups with their corresponding weights for the second orientation
orientation = 2;
for j=IDX(orientation,1):IDX(orientation,2),

    reshape( G(:,j), Height, Width )
    reshape( D(:,j), Height, Width )

end


%==========================================================================

% Same groups with weights (W2), i.e., that depend on the cardinal of the groups
% clear params
% params.cardG_weights = true;

% [G,D, IDX] = get_groups( Height, Width, params );
