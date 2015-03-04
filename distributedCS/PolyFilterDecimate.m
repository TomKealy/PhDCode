%% Filter and decimation
% The function filters a given signal (on a t_axis) using the specified FIR
% filter. The time_axis is truncated by the filter length .
% The output is decimated by factor. The samples and the new
% time_axis are provided.
function [Samples, time_axis] = PolyFilterDecimate(Sig,t_axis,PolyFilter,LengthFilter);
len = length(Sig);
[M, polylen] = size(PolyFilter);

Sig2 = [zeros(1,M-1) Sig];
len2 = length(Sig2);
cosetlen = ceil(len2/M);
x2_pad = [Sig2 zeros(1,cosetlen*M-len2)];
x2_cosets = flipud(reshape(x2_pad,M,cosetlen));
y2 = zeros(1,cosetlen);
for i=1:M
    y2 = y2 + filter(PolyFilter(i,:),1,x2_cosets(i,:));
end
Samples = y2( (polylen+1):(cosetlen-1) );
ind = ((polylen*M+1):M:len) - floor(LengthFilter/2);
time_axis = t_axis(ind);