%% Permutations in Matgraph
% The Matgraph package includes a permutation data type. Here we illustrate
% some of the features of that type.
%% Defining permutations
% Permutations are created using the |permutation| function. The argument
% can either be a positive integer (in which case we create the identity
% permutation or a permuted version of |1:n|. 
p = permutation(6)
q = permutation([1 5 2 4 6 3])
%% Random permutations
% Applying the function |random| to a permutation scrambles its elements.
p = random(permutation(6))
%% Applying a permutation to an element
% If |p| is a permutation and |k| is an integer, then |p(k)| is the result
% of applying |p| to |k|.
disp(q(2))  % This should give 5 since q = (1)(2,5,6,3)(4)
%% Composition
% Composition of permutations is given by the |*| operator.
p = permutation([2 3 4 6 1 5])
q
p*q
q*p
%% Inverses
% The inverse of a permutation can be calculated using |inv|. 
inv(p)
p * inv(p)
%% Powers
% If |p| is a permutation and |k| is an integer, then |p^k| gives the
% |k|-fold composition of |p| with itself. |k| may be negative, giving a
% power of the inverse of |p|. So |p^-1| is the same as |inv(p)|.
p^3
p*p*p
p^-1
inv(p)
%% Conversion to an array
% The |array| method converts a permutation back into an array. 
disp(array(p));
%% Conversion into a list of cycles
% The |cycles| method creates a MATLAB cell array containing the cycles of
% the permutation.
r = random(permutation(17))
cyc = cycles(r);
for k=1:length(cyc)
	disp(['Cycle #', int2str(k), ' is ', int2str(cyc{k})]);
end
%% Conversion into a permutation matrix
% The |matrix| function converts a permutation into a permutation matrix.
q
matrix(q)
%% See also
% Here is a list of all the |permutation| methods.
methods permutation
