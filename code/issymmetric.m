% Checks whether a matrix is symmetric (has to be square)
% Check whether mat=mat^T
% INPUTS: adjacency matrix
% OUTPUTS: boolean variable, {0,1}
% GB, October 1, 2009

function S = issymmetric(mat)

S = false; % default
if mat == transpose(mat); S = true; end