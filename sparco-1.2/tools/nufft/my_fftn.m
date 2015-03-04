function X = my_fftn(x)

if nargin < 1, fftn_test, return, end

X = fft(fft(x).').';
return

X = fft(fft(x, [], 1), [], 2);
return

dd = ndims(x);
X = fft(x, [], 1);
for id=2:dd
	X = fft(X, [], id);
end
return

X = fft(x);
for id=2:dd
	ord = 1:dd;
	ord([1 id]) = [id 1];
	x = permute(X, ord);
	X = fft(x);
	X = ipermute(X, ord);
end


function fftn_test

N = 2^11;
x = randn([N N]);

nrep = 4;
tic
for ii=1:nrep, x = fft(fft(x).').'; end 
printf('inline time %g', toc)

tic
for ii=1:nrep, x = fftn(x); end 
printf('fftn time %g', toc)

tic
for ii=1:nrep, x = my_fftn(x); end 
printf('my_fftn time %g', toc)
