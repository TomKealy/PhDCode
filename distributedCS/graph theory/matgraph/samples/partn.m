%% Partitions in Matgraph
% Many Matgraph functions work with partitions. For example, the |color|
% function returns a partition of the vertex set and |bipmatch| takes a
% partition as one of its arguments. Matgraph provides the |partition| type
% for working with partitions. In all cases, the partition is of the set
% [n]={1,2,...,n} for some n.
%% Creating partitions
% The |partition| function is used to create a new partition.
% |partition(n)| creates a partition of the set [n] in which every element
% is in a part by itself. Alternatively, one can apply |partition| to a
% cell array in which the cells in the array contain the elements in each 
% part.
p = partition(5)
q = partition({[1 2 3],[4 6],[5]})
%% Finding the part that contains a given element
% If |p| is a partition and |k| is an integer, then |p(k)| returns the
% elements in the part of |p| that contains |k|.
disp(q(1))
%% Merging parts
% If |p| is a partition and |j| and |k| are integers, then |merge(p,j,k)|
% is a new partition formed from |p| by combining the parts that contain
% |j| and |k|.
q
merge(q,1,4)
%% Meet and join
% If |p| and |q| are partitions, then |p*q| and |p+q| are the meet and join
% of |p| and |q|, respectively.
p = partition({ [1 3 5], [2 4], [8 7], [6]})
q = partition({ [1 3 2], [4 5], [6 8], [7]})
p*q
p+q
%% Extracting the parts
% |parts(p)| returns a cell array in which each cell contains the elements
% of a part of |p|.
pts = parts(p);
for k=1:length(pts)
	disp(['Part #', int2str(k), ' is ', int2str(pts{k})])
end
%% See also
% Here is a list of all the |partition| methods.
methods partition