function selective(g,n,n0,d)
% selective(g,n,n0,d) --- selective attachment random graph
% overwrite g with a random graph with a degree-d selective attachment
% graph on n vertices starting with n0 isolated vertices.
%
% That is, we begin with n0 isolated vertices. We then add vertices one at
% a time. Each new vertex has d edges to previous vertices in the graph.
% These d edges are drawn at random (without replacement) with the
% probability of joining to a previous vertex u proportion to d(u)+1. The
% process ends when we have n vertices.

clear_edges(g);
resize(g,n);
rmxy(g);

for v=n0+1:n
	wt = deg(g);
	wt = wt + [ones(1,v-1),zeros(1,n-v+1)];
	for k=1:d
		u = weighted_choice( wt);
		add(g,v,u);
		wt(u) = 0;  % to avoid repetition
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function c = weighted_choice(w)
% weighted_choice(w) --- choose an integer at random according to weights.
% w is a list of nonnegative real values
% return an integer c between 1 and length(w) with probability proportional
% to w(c).

n = length(w);
p = cumsum(w);
p = p/p(n);
c = sum(p<rand)+1;
