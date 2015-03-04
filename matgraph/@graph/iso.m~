function [yn,p] = iso(g,h,options)
% [yn,p] = iso(g,h,options) --- is g isomorphic to h?
% Given graphs g and h, determine whether or not the graphs are isomorphic, 
% and if so, return a permutation p such that renumber(g,p) makes g==h.
% 
% Returns yn = 1 if isomorphic and yn = 0 if not. The optional p is a
% permutation giving the isomorphism.
%
% The optional third argument, options, is used to set various parameters
% guiding the search for an isomorphism. Here are the fields that can be
% specified in options and their default values.
%
% options.verbose --- turn off (0, default) or on (1) verbose output.
%
% options.eig --- cospectrality tests. Set this to...
%     negative value: skip this test
%     zero: find characteristic polynomial 
%     positive: calculate eigenvalues and use this value for a threshold
%     for equality.   default value = nv * 20 * eps
%
% options.pretest --- do basic tests such as number of vertices, edges.
%     Set to 1 for on (default) and 0 for off (but don't do that!); these
%     are generally very fast
% 
% options.components --- check that the number and size of components are
%     the same. 1 = check, 0 = don't check. Default is 0
%
% options.distance --- use the distance matrix to distinguish vertex types.
%     1 = check (default), 0 = don't check. Default is 1. This can speed up
%     isomorphism testing for regular graphs.
%

%% Check the options structure and set default values

if nargin < 3
	options.verbose = 0;
end

if ~isfield(options,'verbose')
	options.verbose = 0;
end

if ~isfield(options,'eig')
	options.eig = nv(g) * 20 *eps;
end

if ~isfield(options,'pretest')
	options.pretest = 1;
end

if ~isfield(options,'components')
	options.components = 0;
end

if ~isfield(options,'distance')
	options.distance = 1;
end

if ~isfield(options,'debug')
	options.debug = false;
end


verb = options.verbose;

if verb
	disp('Searching for an isomorphism')
end

% default in case the graphs are not isomorphic
yn = false;
p = permutation(0);

%%
% We begin with a set of basic checks that can easily determine if the
% graphs are not isomorphic.


% check if the graphs are simply equal

if verb
	disp('Checking if the graphs are identical')
	tic
end

if (g==h)
	yn = true;
    p = permutation(nv(g));
	if verb
		toc
	end
    return
end


% number of vertices & edges
if options.pretest
	if ~iso_pretest(g,h,options)
		if verb
			toc
		end
		return
	end
end


%% Create vertex distinguishing information. This is kept in n-row matrices
% Mg and Mh.

% first set of distinguishing mark is the degree sequence of its neighbors

if verb
	disp('Calculating vertex neighbor degree sequences')
end

n = nv(g);

Mg = deg(g);
Mg = Mg(:);
Mh = deg(h);
Mh = Mh(:);
	
Mg = refine_types(g,Mg);
Mh = refine_types(h,Mh);

% check if they're the same

if verb
	toc
end

if any(any( sortrows(Mg) ~= sortrows(Mh) ))
	return
end

% the second set of distinguishing marks is from the distance matrix

if options.distance
	if verb
		disp('Calculating distance matrices')
	end
	Dg = sort(dist(g),2);
	Dh = sort(dist(h),2);

	% append to Mg and Mh

	Mg = [Mg, Dg(:,2:end)];
	Mh = [Mh, Dh(:,2:end)];
	
	if verb
		toc
	end
	
end


% make sure they're the same
if any(any( sortrows(Mg) ~= sortrows(Mh) ))
	if options.verbose
		disp('Graphs have different distance structure')
	end
	return
end


%% Now we assign a "type" to each vertex of the graph based on its row in
% Mg (or Mh). Vertices of different types can NOT be linked by an
% isomorphism.


% determine type number for each vertex

if verb
	disp('Refining')
end

gtypes = matrix2types(Mg);
htypes = matrix2types(Mh);

gtypes = ultra_refine_types(g,gtypes);
htypes = ultra_refine_types(h,htypes);


ntypes = max(gtypes);

if any(sort(gtypes) ~= sort(htypes))
	yn = 0;
	if verb 
		toc
	end
	return
end


if verb
	disp(['Number of different types of vertices = ', int2str(ntypes)])
	toc
end

[new_gtypes,new_htypes] = iso_match(g,h,gtypes,htypes,options);

if max(new_gtypes) < n || max(new_htypes) < n
	return
end

p1 = permutation(new_gtypes);
p2 = permutation(new_htypes);

p = inv(p2)*p1;

gcopy = graph;
copy(gcopy, g);
renumber(gcopy,p);
if gcopy==h
	yn = 1;
	if options.verbose
		disp('Graphs are isomorphic');
		toc
	end
else
	yn = 0;
	p = permutation(0);
	if options.verbose
		disp('Graphs are NOT isomorphic');
		toc
	end
end
free(gcopy);


%% main work engine
function [new_gtypes, new_htypes] = iso_match(g,h,gtypes,htypes,options)

ntypes = max(gtypes); 
n = nv(g);

if ntypes == n
	new_gtypes = gtypes;
	new_htypes = htypes;
	return
end

% find a g-vertex of the most frequent type

type_count = zeros(1,ntypes);

% find a vertex of the most numerous type
for k=1:ntypes
	type_count(k) = length(find(k==gtypes));
end
max_idx = find(type_count == max(type_count));
max_type = max_idx(1);

% get a g-vertex of type max_type; call it v
	
v = find(gtypes == max_type);
v = v(1);

new_num = ntypes + 1; % create a new type for v
gtypes(v) = new_num;
new_gtypes = ultra_refine_types(g,gtypes);


% find all h vertices of the same type as v
% call that list wlist

wlist = find(htypes == max_type);

nw = length(wlist);
for i = 1:nw
	% check if G-v is plausible match to H-w
	w = wlist(i);
	
	if options.verbose
		disp(['Trying ', int2str(v), ' --> ', int2str(w)])
		toc
	end
	
	if options.pretest
		g_v = graph; copy(g_v,g); delete(g_v,v);
		h_w = graph; copy(h_w,h); delete(h_w,w);
		if (~iso_pretest(g_v,h_w,options))
			if options.verbose
				disp('failed pretest')
			end
			free(g_v);
			free(h_w);
			continue		
		end
		free(g_v);
		free(h_w);
	
		if options.verbose
			disp('Pretest passed')
		end
	end

	% plausible match, so assign new type to w and recurse
	htypes_tmp = htypes;
	htypes_tmp(w) = new_num;
	htypes_tmp = ultra_refine_types(h,htypes_tmp);
	
	
	% htypes_tmp should be same as gtypes (up to sort order)
	if any(sort(htypes_tmp) ~= sort(new_gtypes))
		if options.verbose
			disp('Attempted match failed')
		end
		continue;
	end
	
	
	if options.verbose
		disp(['Number of vertex types is now ', int2str(max(htypes_tmp))])
		toc
	end
	
	% recurse
	[new_gtypes, new_htypes] = iso_match(g,h,new_gtypes,htypes_tmp,options);
	
	if (max(new_gtypes) == n) && (max(new_htypes)==n)
		p1 = permutation(new_gtypes);
		p2 = permutation(new_htypes);

		p = inv(p2)*p1;

		gcopy = graph;
		copy(gcopy, g);
		renumber(gcopy,p);
		if gcopy==h
			yn = 1;
		else
			yn = 0;
		end
		free(gcopy);
		
		if yn
			return
		end
	end
end

new_gtypes = gtypes;
new_htypes = htypes;
	

%% matrix2types -- convert rows of M to type numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tlist = matrix2types(M)
% convert the rows of M into a column of row types

T = unique(sortrows(M),'rows');

n = size(M,1);

% determine type number for each row

tlist = zeros(n,1);

for v=1:n
	r = M(v,:);
	tlist(v) = find_row(r,T);
end



%% refine_types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newtypes = refine_types(g,types)
% newtypes = refine_types(g,types) -- refine types based on neighbors' types
% given an n-vector of vertex types, augment this by considering the types
% of the neighbors


maxd = max(deg(g));
n = nv(g);

M = zeros(n,maxd+1);

for v=1:n
	Nv = neighbors(g,v);
	row = types(Nv);
	row = sort(row);
	row = row(:)';  % make sure it's a row vector
	row = [types(v), row, zeros(1, maxd-length(row))];
	M(v,:) = row;
end


newtypes = matrix2types(M);

%% ultra_refine_types -- repeatedly apply refine_types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newtypes = ultra_refine_types(g,types)
% repeatedly apply refine_types until no change

while true
	newtypes = refine_types(g,types);
	if all(newtypes == types)
		return
	end
	types = newtypes;
end



%% iso_pretest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yn = iso_pretest(g,h,options)
% basic tests that can rule out isomorphism

% # of vertices and edges

yn = 0;
if any(size(g) ~= size(h))
	return
end

% degree sequences

dg = sort(deg(g));
dh = sort(deg(h));
if any(dg ~= dh)
	if options.verbose
		disp('Graphs have different degree sequences')
	end
	return
end  

% components

if options.components


	gc = components(g);
	hc = components(h);

	if np(gc) ~= np(hc)
		return
	end

	%  components must all have the same sizes
	gcp = parts(gc);
	hcp = parts(hc);

	gx = zeros(1,np(gc));
	gy = gx;

	for k=1:np(gc)
		gx(k) = length(gcp{k});
		gy(k) = length(hcp{k});
	end
	gx = sort(gx);
	gy = sort(gy);
	if any(gx ~= gy)
		if options.verbose
			disp('Graphs have different component sizes')
		end
		return
	end
end % end component subtest

% cospectrality tests

if options.eig == 0
	
	A = double(matrix(g));
	B = double(matrix(h));
	pA = round(poly(A));
	pB = round(poly(B));
	if any(pA ~= pB)
		if options.verbose
			disp('Graphs have different characteristic polynomials')
		end
		return;
	end
end


if options.eig > 0
	eA = sort(eig(g));
	eB = sort(eig(h));
	diff = norm(eA-eB);
	if (diff > options.eig)
		if options.verbose
			disp(['Graphs have different spectra (tol = ', ... 
				num2str(options.eig),')']);
		end
		return
	end
end


yn = 1;
