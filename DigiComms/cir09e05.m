%dc09e05.m: A Linear Block Code in Example 9.5
% Constructs the codeword matrix and finds its minimum distance.
clear
K=4;  L=2^K; % Message size and the number of codewords
for i=1:L,  M(i,:)=deci2bin1(i-1,K);  end
M % A series of K-bit binary numbers
% Generator matrix
G=[1 1 0 1 0 0 0; 0 1 1 0 1 0 0; 1 1 1 0 0 1 0; 1 0 1 0 0 0 1];
% To generate codewords
Codewords =rem(M*G,2) % Modulo-2 multiplication Eq.(9.4.5)
% Find the minimum distance by Eq.(9.4.8)
Minimum_distance =min(sum((Codewords(2:L, :))'))