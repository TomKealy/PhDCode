function r = regop(y,p,adj,modul,bc)

if ~exist('bc','var')
	bc = 'zero';
end;

if ~exist('adj','var')
    adj = false;
end;

if ~exist('modul','var')
    modul = false;
end;

% XXX
if ~modul
    fdfwd = poly(exp(p(p<=0)));
    fdbck = poly(exp(p(p>0)));
else
    fdfwd = poly(-exp(p(p<=0)));
    fdbck = poly(-exp(p(p>0)));
end;

r = y;

if ~adj
        r = fdfwd(1)*r;
	for ii=2:length(fdfwd)
    	r = r + fdfwd(ii) * shift(r,ii-1,bc);
	end;
	
        r = fdbck(1)*r;
	for ii=2:length(fdbck)
		r = r + fdbck(ii) * shift(r,-ii+1,bc);
	end;
else
        r = fdfwd(1)*r;
	for ii=2:length(fdfwd)
		r = r + conj(fdfwd(ii)) * shift(r,-ii+1,bc);
	end;

        r = fdbck(1)*r;
	for ii=2:length(fdbck)
		r = r + conj(fdbck(ii)) * shift(r,ii-1,bc);
	end;	
end; 

