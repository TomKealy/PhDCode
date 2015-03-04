function data = prob055(varargin)
%PROB055  Mondrian.
%
%   PRELIMINARY! DON'T USE.
%
%   Examples:
%   P = prob055;  % Creates the default 055 problem.
%
%   References:
%
%   [TsaiDono:2006] Y. Tsaig and D.L. Donoho, Extensions of
%     compressed sensing, Signal Process., 86:3 (2006), pp. 549-571.
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: probMondrian.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{});

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',0); end;

% Set up the data
n      = 1024;
blocks = [ 75, 32,414,271,0.54424794766829; ...
          583,261,744,338,0.86735677370723; ...
          606,502,873,907,0.32314338345102; ...
          852,161,909,314,0.98745160741150; ...
          296,726,724,831,0.71725830291318; ...
          534,769,712,918,0.01080709434893];

signal = zeros(n,n);
for i=1:size(blocks,1)
   r0 = blocks(i,1);   r1 = blocks(i,3);
   c0 = blocks(i,2);   c1 = blocks(i,4);
   signal(r0:r1,c0:c1) = signal(r0:r1,c0:c1) + blocks(i,5);
end

% Set up the problem
data.signal   = signal;
data.op.Dirac = opDirac(1);
data.B        = data.op.Dirac;
data.b        = 1;
data          = completeOps(data);

% Additional information
info.title           = 'Mondrian image';
info.thumb           = 'figProblem055';
info.citations       = {'TsaiDono:2006'};
info.fig{1}.title    = 'Original image';
info.fig{1}.filename = 'figProblem055Mondrian';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  image(data.signal * 59.5), colormap gray;
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)
  
  if opts.update
     idx = round(linspace(1,n+1,129));
     idx = idx(1:end-1);
     P   = data.signal(idx,idx);
     P   = (P(1:2:end,:) + P(2:2:end,:)) / 2;
     P   = (P(:,1:2:end) + P(:,2:2:end)) / 2;
     thumbwrite(P, info.thumb, opts);
  end
end
