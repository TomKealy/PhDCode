function indicator = is_contained(A,B)
% test if A is contained in B
ContainedNum = length(intersect(A,B));
indicator = (ContainedNum == length(A));