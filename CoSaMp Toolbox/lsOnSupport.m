function alpha = lsOnSupport(y,PhiSub,normBound)
%% Given measurements y and submatrix PhiSub, recover alpha
%
% Usage: alpha = lsOnSupport(y,PhiSub,normbound)
% Inputs: y - measurements obtained
%         PhiSub - submatrix of Phi
%         normBound - an upper bound on the norm of the alpha to be
%                     recovered for use in Tikhonov regularization
%
% Outputs: alpha - recovered coefficients
%
% Most recent change - 7/31/2012
%
% Copyright 2012, Mark Davenport, Deanna Needell, Michael Wakin
%
% This file is part of Signal Space CoSaMP Toolbox version 1.0.
%
%    Signal Space CoSaMP Toolbox is free software: you can redistribute it 
%    and/or modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, either version 3 of the 
%    License, or (at your option) any later version.
%
%    Signal Space CoSaMP Toolbox is distributed in the hope that it will be 
%    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Signal Space CoSaMP Toolbox.  If not, 
%    see <http://www.gnu.org/licenses/>.

%% Recover alpha via Tikhonov regularized least-squares
[U,s,V] = csvd(PhiSub);
[alpha2,lambda,maxIterReached] = lsqi(U,s,V,y,normBound);
%alpha = pinv(PhiSub)*y;
%alpha = PhiSub\y;

% We have modified lsqi.m to return a flag if it does not converge (rather
% than printing a warning). 
if maxIterReached
    % lsqi did not converge; use cvx (slower) instead
    [aa,bb] = size(PhiSub);
    cvx_begin quiet
      variable alpha(bb) complex;
      minimize norm(y - PhiSub*alpha)
      subject to
      norm(alpha) <= normBound;
    cvx_end
    
    % disp(['normBound = ' num2str(normBound)]);
    % disp([' norm(alpha) = ' num2str(norm(alpha)) ', norm(y - PhiSub*alpha) = ' num2str(norm(y-PhiSub*alpha))]);
    % disp([' norm(alpha2) = ' num2str(norm(alpha2)) ', norm(y - PhiSub*alpha2) = ' num2str(norm(y-PhiSub*alpha2))]);
    % disp([' norm(alpha-alpha2) = ' num2str(norm(alpha-alpha2))]);
else   
    alpha = alpha2;
end

