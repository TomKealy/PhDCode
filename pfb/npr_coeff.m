%NPR_COEFF generates near NPR filter bank coefficients
%  COEFF = NPR_COEFF(N,L) generates the filter coefficients
%  for a near perfect reconstruction filter bank.
%  The number of subbands will be N, and L is the number of 
%  filter taps used per subband. The output COEFF will have 
%  size (N/2,L).
%
%  The prototype is constructed starting with an equiripple
%  approximation to a 'root raised error function' shaped filter.
%
%  NPR_COEFF(N,L) with no output arguments plots the magnitude
%  and the prototype filter in the current figure window.
%
%  See also npr_analysis, npr_synthesis, npr_coeff
%
% (c) 2007 Wessel Lubberhuizen
function coeff = npr_coeff(N,L,K)

if ~exist('N','var')
    N=256;  % if the number of subband is not specified, use a default value
end

if ~exist('L','var');
    L=128; % if the number of taps is not specified, use a default value
end

if ~exist('K','var');
    % if the K value is not specified us a default value.
    % these values minimize the reconstruction error.
    switch(L)
        case 8,    K=4.853;
        case 10,   K=4.775;
        case 12,   K=5.257;
        case 14,   K=5.736;
        case 16,   K=5.856;
        case 18,   K=7.037;
        case 20,   K=6.499;
        case 22,   K=6.483;
        case 24,   K=7.410;
        case 26,   K=7.022;
        case 28,   K=7.097;
        case 30,   K=7.755;
        case 32,   K=7.452;
        case 48,   K=8.522;
        case 64,   K=9.396;
        case 96,   K=10.785;
        case 128,  K=11.5; 
        case 192,  K=11.5;
        case 256,  K=11.5;
        otherwise, K=8;
    end
end

% divide the number of subbands by two, because we're creating two
% overlapping filterbanks.
M = N /2;

F= (0:(L*M-1))/(L*M);


% The prototype is based on a root raised error function
A = rrerf(F,K,M);
N=length(A);

n=0:(N/2-1);
A(N-n)=conj(A(2+n));
A(1+N/2)=0;

B=ifft(A);
B=fftshift(B);
A=fftshift(A);

if nargout==0
    % if not output arguments are specified, 
    % create some plots to visualize the result.
      
    figure(1);
    subplot(2,1,1);
    plot(20*log10(abs(A)));
    title('prototype filter - frequency response')
    
    grid on;
    axis([1 length(A) -350 0]);
    subplot(2,1,2);
    plot(20*log10(abs(B)))
    title('prototype filter - impulse response')
    ylabel('dB')
    grid on;
    axis([1 length(B) -350 0]);

    N=16384;
    F=(0:N-1)/N;
    A=rrerf(F,K,M);
    N=length(A);
    n=0:(N/2-1);
    A(N-n)=conj(A(2+n));
    A(1+N/2)=0;
    A=fftshift(A);

    F=F-0.5;
    W=2*pi*F;
    H=freqz(B,1,W);
    figure(2);
    subplot(2,1,1);
    plot(2*F*M,20*log10(abs([H' A'])));
    axis([-M M -400 0]);
    xlabel('frequency (channel)');
    ylabel('power (dB)');
    grid on;

    subplot(2,1,2);
    plot(2*F*M,20*log10(abs(abs(H)-abs(A))));
%    xlabel('frequency (channel)');
    ylabel('deviation (dB)');
    grid on;

    H1=H(      1:N-N/M/1);
    H2=H(1+N/M/2:N-N/M/2);
    H3=H(1+N/M/1:N      );   
    
    

    F=F(1:N-N/M/1);
    subplot(2,1,2);
    plot(2*F*M+1,20*log10(abs([H1.*H1;H2.*H2;H3.*H3;H1.*H2;H2.*H3;H1.*H3])));
    xlabel('frequency (channel)');
    ylabel('reconstr. error (dB)');

    grid on;
    coeff=[];
else
    B=B/sum(B);
    coeff = reshape(B,M,L);
end


function A = rrerf(F,K,M)

x = K*(2*M*F-0.5);
A= sqrt(0.5*erfc(x));