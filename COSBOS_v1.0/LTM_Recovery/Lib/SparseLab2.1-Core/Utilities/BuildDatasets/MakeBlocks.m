function sig = MakeBlocks(n)

t = (1:n) ./n;
pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
sig = zeros(size(t));
for j=1:length(pos)
    sig = sig + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
end
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
