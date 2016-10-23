%NPR_ANALYSIS Near perfect reconstruction analysis
%  Y = NPR_ANALYSIS(COEFF,X) separates the input signal X into
%  subbands. Each subband is a row in X. COEFF is a two dimensional array, 
%  containing the filter coefficients. The number of rows in X will be
%  be twice the number of rows in COEFF. 
%
%  See also npr_synthesis, npr_coeff
%
% (c) 2007 Wessel Lubberhuizen


function y = npr_analysis(coeff,x)

% number of channels
N=size(coeff,1);

% number of slices
M=ceil(length(x)/N);

% create polyphase input signals
x1=reshape(x,N,M);
x2=x1;
for i=1:N;
    x2(i,:) = x2(i,:) * exp(sqrt(-1)*pi*(i-1)/N);
    x2(i,2:2:M) = -x2(i,2:2:M);
end
% apply channel filters
coeff = fliplr(coeff);
for i=1:N
    x1(i,:) = filter(coeff(i,:),1,x1(i,:));
    x2(i,:) = filter(coeff(i,:),1,x2(i,:));    
end

% apply dft
x1 = ifft(x1,[],1)*N;
x2 = ifft(x2,[],1)*N;

% assemble even and odd channels
y = [reshape(x1,1,N,M); reshape(x2,1,N,M)];
y = reshape(y,2*N,M);
