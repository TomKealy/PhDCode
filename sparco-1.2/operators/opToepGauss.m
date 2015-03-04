function [op,d] = opToepGauss(m,n,type,normalized)
%OPTOEPGAUSS  Toeplitz matrix with Gaussian entries
%
%   OP = OPTOEPGAUSS(M,N,TYPE,NORMALIZED) creates an M by N
%   Toeplitz matrix with Gaussian entries. TYPE can either be
%   'toeplitz' or 'circular'. For the 'toeplitz' type matrix, m+n-1
%   different generating entries are generated at random, whereas
%   for the 'circular' matrix, only max(m,n) are needed. When the
%   TYPE field is empty [], or not specified 'toeplitz' is chosen
%   by default. Setting the NORMALIZED flag scales the columns of
%   the Toeplitz matrix to unit norm. Multiplication is implemented
%   using the fast Fourier transform. The use of such matrices in
%   compressed sensing was suggested by:
%
%   [1] W. U. Bajwa, J. D. Haupt, G. M. Raz, S. J. Wright, and
%       R. D. Nowak, Toeplitz-Structured Compressed Sensing Matrices,
%       IEEE/SP 14th Workshop on Statistical Signal Processing,
%       pp. 294-298, August 2007.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opToepGauss.m 1052 2008-07-02 20:52:08Z ewout78 $

if nargin < 3 || isempty(type)
    type = 'Toeplitz';
end

if nargin < 4
   normalized = 0;
end


switch lower(type)
    case 'circular'
       % Generate the entries of the matrix
       k  = max(m,n);
       d  = randn(k,1);
       df = fft(d);

       if normalized
          if m == n
             s = 1 / norm(d);
          else
             s = zeros(n,1);
             for i=1:n
               idx = [1:m-i+1,k-i+2:k-max(0,i-m-1)];
               s(i) = 1 / sqrt(sum(d(idx).^2));
             end             
          end         
       else
         s = 1;
       end
       
       op = @(x,mode) opToepGaussCircular_intrnl(df,s,m,n,x,mode);
       
    case 'toeplitz'
       % Generate the entries of the matrix
       k  = max(m,n);
       d  = randn(2*k,1); d(m+1:end-n+1) = 0;
       df = fft(d);

       if normalized
          s = zeros(n,1); k = length(d);
          for i=1:n
             idx = [1:m-i+1,k-i+2:k-max(0,i-m-1)];
             s(i) = 1 / sqrt(sum(d(idx).^2));
          end             
       else
         s = 1;
       end

       op = @(x,mode) opToepGauss_intrnl(df,s,m,n,x,mode);
    otherwise
        error('Unrecognized type parameter');
end


function y = opToepGauss_intrnl(df,s,m,n,x,mode)
checkDimensions(m,n,x,mode);
if mode == 0
    y = {m,n,[0,1,0,1],{'ToepGauss','toeplitz'}};
elseif mode == 1
    k = max(m,n);
    y = opToepGaussCircular_intrnl(df,1,2*k,2*k,[s.*x;zeros(2*k-n,1)],mode);
    y = y(1:m);
else
    k = max(m,n);
    y = opToepGaussCircular_intrnl(df,1,2*k,2*k,[x;zeros(2*k-m,1)],mode);
    y = s.*y(1:n);  
end


function y = opToepGaussCircular_intrnl(df,s,m,n,x,mode)
checkDimensions(m,n,x,mode);
if mode == 0
    y = {m,n,[0,1,0,1],{'ToepGauss','circular'}};
elseif mode == 1
    y = ifft(df.*fft([s.*x;zeros(max(m,n)-n,1)]));
    y = y(1:m);
    if isreal(x), y = real(y); end;
else
    y = ifft(conj(df).*fft([x;zeros(max(m,n)-m,1)]));
    y = s.*y(1:n);
    if isreal(x), y = real(y); end;
end
