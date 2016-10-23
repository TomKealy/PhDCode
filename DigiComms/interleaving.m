function y=interleaving(x,Nrow,Ncbps,Nbpsc,Ksee)
% x    : Frame of bits/symbols to be interleaved
% Nrow : Number of rows of the block interleaver
% Ncbps: Number of coded bits per symbol
% Nbpsc: Number of coded bits per subcarrier
% y    : Interleaved frame of bits/symbols
y=x;   % To make y have the same size as x
k=0:Ncbps-1;
i = Ncbps/Nrow*mod(k,Nrow) + floor(k/Nrow); % 1st permutation
% Ncol=Ncbps/Nrow; i=Nrow*mod(k,Ncol) + floor(k/Ncol); % 1st permutation
z=x(i+1); i=k; s= max(Nbpsc/2,1); 
if nargin>4&Ksee>0,  Ncol=Ncbps/Nrow; Z=reshape(z,Nrow,Ncol),  end 
j=s*floor(i/s) + mod(i+Ncbps-floor(Nrow*i/Ncbps),s); % 2nd permutation
y=z(j+1);