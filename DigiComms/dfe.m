function  [c,b]=dfe(c,b,ys,ds,e,delta)
% updates decision feedback equalizer(DFE) coefficients
%Input:  ys= History of channel output 
%        ds= Data sequence or detector output history
%        e = Error (discrepancy between the DTR input and output)
%        delta= Step size 
%Output: c,b= Updated DFE coefficients
c=c+delta*e*ys; % Eq.(6.2.17a)
b=b-delta*e*ds; % Eq.(6.2.17b)
