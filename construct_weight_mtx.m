function weights = construct_weight_mtx(weight_list, Adj)

weights = zeros(size(Adj));
positions = find(tril(Adj));

weights(positions) = weight_list;

weights = weights - weights.';   %// ' is complex conjugate; not a big deal here, but something to know

find(Adj) == find(weights);   %// Not sure what this is meant to do; maybe an assert?

end 