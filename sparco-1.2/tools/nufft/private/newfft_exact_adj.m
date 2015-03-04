 function x = newfft_exact_adj(st, X, om)
% forward
nthread = 1; % todo
useloop = false; % todo

if ~isempty(st.om)
	om = st.om;
end
if exist('dtft_mex') == 3
	x = jf_mex('dtft,adjoint', double(om'), double(X), ...
		int32(st.Nd), int32(nthread));
	if any(st.n_shift)
		error 'n_shift not done'
	end
else
	x = dtft_adj(X, om, st.Nd, st.n_shift, useloop);
end
