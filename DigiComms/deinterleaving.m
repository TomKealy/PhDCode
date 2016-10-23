function x=deinterleaving(y,Nrow,Ncbps,Nbpsc)
% y    : Frame of bits/symbols to be deinterleaved
% Nrow : Number of rows of the block interleaver
% Ncbps: Number of coded bits per symbol
% Nbpsc: Number of coded bits per subcarrier
% x    : Deinterleaved frame of bits/symbols 
x=y; % To make x have the same size as y
i=0:Ncbps-1; s= max(Nbpsc/2,1); 
j=s*floor(i/s) + mod(i+Ncbps-floor(Nrow*i/Ncbps),s); % 1st permutation
z(j+1)=y;  k=i;
i = Ncbps/Nrow*mod(k,Nrow) + floor(k/Nrow); % 2nd permutation
%Ncol=Ncbps/Nrow; i=Nrow*mod(k,Ncol) + floor(k/Ncol); % 2nd permutation
x(i+1)=z;