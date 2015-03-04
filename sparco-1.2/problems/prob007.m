function data = prob007(varargin)
%PROB007  Gaussian ensemble, sign/spikes signal.
%
%   PROB007 creates a problem structure.  The generated signal will
%   consist of K = 20 sign spikes and have length N = 2560. The
%   signal is measured using an M = 600 by N Gaussian ensemble with
%   orthogonal rows.
%
%   The following optional arguments are supported:
%
%   PROB007('k',K,'m',M,'n',N,'scale',SCALE,flags) is the same as
%   above, but with a measurement matrix of size SCALE*M by SCALE*N
%   and a signal length of SCALE*N. The number of spikes in the signal
%   is set to SCALE*K. When scale is not an integer, all parameters are
%   rounded to the nearest integer after scaling. By default, the SCALE
%   is set to 1. The 'noseed' flag can be specified to suppress
%   initialization of the random number generators. Each of the
%   parameter pairs and flags can be omitted.
%
%   Examples:
%   P = prob007;  % Creates the default 007 problem.
%
%   References:
%
%   [l1magic] E.J. Candes and J. Romberg, l1-Magic, 2005
%     http://www.l1-magic.org/ 
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob007.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'k','m','n','scale'});
scale       = getOption(parm,'scale', 1);
k           = getOption(parm,'k',    20);
m           = getOption(parm,'m',   600);
n           = getOption(parm,'n',  2560);
info.name   = 'sgnspike';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',0); end;

% Set up the data
k         = max(1,round(k * scale));
m         = max(1,round(m * scale));
n         = max(1,round(n * scale));
p         = randperm(n); p = p(1:k);
signal    = zeros(n,1);
signal(p) = sign(randn(k,1));

% Set up the problem
data.signal      = signal;
data.op.Gaussian = opGaussian(m,n,4);
data.M           = data.op.Gaussian;
data.x0          = data.signal;
data.b           = data.M(data.x0,1);
data             = completeOps(data);

% Additional information
info.title           = 'Sign spikes, real domain';
info.thumb           = 'figProblem007';
info.citations       = {'CandRomb:2007'};
info.fig{1}.title    = 'Sign spike signal';
info.fig{1}.filename = 'figProblem007Signal';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(1:n,data.signal,'b-');
  xlim([1,n]); ylim([-1.2,1.2]);
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)
  
  if opts.update
     P = ones(128,128,3);
     P = thumbPlot(P,1:n,data.signal,[0,0,1]);
     P = (P(1:2:end,:,:) + P(2:2:end,:,:)) / 2;
     P = (P(:,1:2:end,:) + P(:,2:2:end,:)) / 2;
     thumbwrite(P, info.thumb, opts);
  end
end
