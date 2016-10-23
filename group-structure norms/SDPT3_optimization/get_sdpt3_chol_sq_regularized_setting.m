function [ At, b, c, blk ] = get_sdpt3_chol_sq_regularized_setting( XtX_over_n, XtY_over_n, G, D, lambda )

%GET_SDPT3_SETTINGS

[p, ng] = size(G);

group_size     = sum( G, 1 );
sum_group_size = sum( group_size );

% max b^T u
% s.t. c-At u \in Cone


n_row_At       = (2+p)+(ng+sum_group_size)+3;
n_variables    = 1+1+ng+p;
n_nonzero_elts = (2+p*p)+(ng+sum_group_size)+(2+ng);

blk = cell(1, 2);
At  = cell(1, 1);
c   = cell(1, 1);

% form vector b

b = [-1 ; -lambda ; zeros(ng,1) ; XtY_over_n];

% form vector c (WARNING, use rotated cone)
% 1:ng+sum_group_size, lorentz cones for ||w_G|| <= u_G
% ng+sum_group_size+1: ng+sum_group_size+3, rotated cone for Omega^2 <= u_2
% end, rotated cone for ||Rw||^2 <= u_1

c{1} = sparse( [ng+sum_group_size+1, ng+sum_group_size+2, ng+sum_group_size+4, ng+sum_group_size+5], ...
               [1, 1, 1, 1], [1, 1, 1, 1]/sqrt(2), n_row_At, 1, 4);


blk{1,1} = 'q';
blk{1,2} = [1+group_size, 3, 2+p];



ROW_INDEX = zeros(n_nonzero_elts,1);
COL_INDEX = zeros(n_nonzero_elts,1);
VALUE     = zeros(n_nonzero_elts,1);

t=1;
end_idx = ng+3:n_variables;

for g=1:ng,
       
    ROW_INDEX(t)                   = t;
    ROW_INDEX(t+1:t+group_size(g)) = t+1:t+group_size(g);
    
    COL_INDEX(t)                   = g+2;
    COL_INDEX(t+1:t+group_size(g)) = end_idx( G(:,g) );
    
    VALUE(t)                   = -1;
    VALUE(t+1:t+group_size(g)) = -D( G(:,g), g);
    
    t = t + group_size(g) + 1;
    
end

%--------------------------------------------------------------------------

ROW_INDEX(t) = t;
COL_INDEX(t) = 2;
VALUE(t)     = -1/sqrt(2);

t = t+1;

ROW_INDEX(t) = t;
COL_INDEX(t) = 2;
VALUE(t)     = 1/sqrt(2);

t = t+1;

ROW_INDEX(t:t+ng-1) = t;
COL_INDEX(t:t+ng-1) = 3:2+ng;
VALUE(t:t+ng-1)     = -1;

t = t+ng;

%--------------------------------------------------------------------------

ROW_INDEX(t) = t-ng+1;
COL_INDEX(t) = 1;
VALUE(t)     = -1/sqrt(2);

t = t+1;

ROW_INDEX(t) = t-ng+1;
COL_INDEX(t) = 1;
VALUE(t)     = 1/sqrt(2);

%--------------------------------------------------------------------------

[tmp1, tmp2] = ndgrid(n_row_At-p+1:n_row_At, n_variables-p+1:n_variables);

R = chol( XtX_over_n );

ROW_INDEX(t+1:end) = tmp1(:);
COL_INDEX(t+1:end) = tmp2(:);
VALUE(t+1:end)     = -R(:);

At{1} = sparse( ROW_INDEX, COL_INDEX, VALUE, n_row_At, n_variables, n_nonzero_elts);

end
