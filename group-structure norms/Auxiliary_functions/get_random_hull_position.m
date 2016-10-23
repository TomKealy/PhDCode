function new_hull_pattern = get_random_hull_position( hull_pattern, height, width )

%GET_RANDOM_HULL_POSITION
% hull_pattern is the pattern we want to find a random position for
% By default, hull_pattern is located in the top left corner
%
%   Copyright (c) 2009 Rodolphe Jenatton. All rights reserved.

[row_idx,col_idx] = ind2sub([height, width], hull_pattern);

row_translation = floor( ( height - max(row_idx) + 1)*rand() );
col_translation = floor( ( width -  max(col_idx) + 1)*rand() );

row_idx = row_translation + row_idx;
col_idx = col_translation + col_idx;

new_hull_pattern = false( height * width,1);
new_hull_pattern( sub2ind( [height, width], row_idx, col_idx) ) = true;


end
