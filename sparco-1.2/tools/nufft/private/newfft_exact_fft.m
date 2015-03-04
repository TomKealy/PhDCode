 function X = newfft_exact_fft(st, x, om)
% forward
nthread = 1; % todo
useloop = false; % todo

if ~isempty(st.om)
	om = st.om;
end
if exist('dtft_mex') == 3
	X = jf_mex('dtft,forward', double(om'), double(x), int32(nthread));
	if any(st.n_shift)
		error 'n_shift not done'
	end
else
	X = dtft(x, om, st.n_shift, useloop);
end
