%% Test sense.m
% 
% Most recent change - 9/9/2011
%
% Copyright 2011, Mark Davenport, Michael Wakin
%
% This file is part of DPSS Approximation and Recovery Toolbox version 1.0.
%
%    DART is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    DART is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DART.  If not, see <http://www.gnu.org/licenses/>.

%% Initialize
clear all; close all;

addpath('../')

N = 12;
M = 3;
K = 2;

x = zeros(N,1); % Generate test signal
x(1:K) = randn(K,1);

%% Test random demodulator
senseOpts.sampler = 'rdmdemod';  % Select random demodulator
senseOpts.pnseq = sign(randn(N,1)); % Choose random pn sequence
senseOpts.subsamp = N/M; % Subsample by factor of N/M

% Generate Phi
pnseq = senseOpts.pnseq;
subsamp = senseOpts.subsamp;
Phi = zeros(M,N);
for jj=1:M,
    Phi(jj,(jj-1)*subsamp+1:jj*subsamp) = pnseq((jj-1)*subsamp+1:jj*subsamp);
end

% Compare sense with Phi
y1 = Phi*x;
y2 = sense(x,0,senseOpts);
x1 = Phi'*y1;
x2 = sense(y1,1,senseOpts);
disp(['Random demodulator: ' num2str(norm(y1-y2)) ' - ' num2str(norm(x1-x2))])

%% Test random sampler
senseOpts.sampler = 'rdmsample';  % Select random demodulator
senseOpts.idx = [1:M];

% Generate Phi
Phi = zeros(M,N);
for jj=1:M,    
    Phi(jj,senseOpts.idx(jj)) = 1;
end

% Compare sense with Phi
y1 = Phi*x;
y2 = sense(x,0,senseOpts);
x1 = Phi'*y1;
x2 = sense(y1,1,senseOpts);
disp(['Random sampler: ' num2str(norm(y1-y2)) ' - ' num2str(norm(x1-x2))])