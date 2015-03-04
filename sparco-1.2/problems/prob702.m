function data = prob702(varargin)
%PROB702  GPSR example: Dirac basis, 2-dimensional sign/spikes signal.
%
%   PROB702 creates a problem structure.  The generated signal will
%   consist of N = 64 by N image with +1/-1 valued impulse noise at
%   each pixel occurring with probability P = 0.01. The signal is
%   blurred by convolution with an 8 by 8 blurring mask and normally
%   distributed noise with standard deviation SIGMA = 0.01 is added
%   to the final signal.
%
%   The following optional arguments are supported:
%
%   PROB702('n',N,'p',P,'sigma',SIGMA,flags) is the same as above, but
%   with signal size N by N, spikes occurring with probability P and
%   an additive noise level of SIGMA. The 'noseed' flag can be specified
%   to suppress initialization of the random number generators. Each of
%   the parameter pairs and flags can be omitted.
%
%   Examples:
%   P = prob702;  % Creates the default 702 problem.
%
%   References:
%
%   [FiguNowaWrig:2007] M. Figueiredo, R. Nowak and S.J. Wright,
%     Gradient projection for sparse reconstruction: Application to
%     compressed sensing and other inverse problems, Submitted,
%     2007. See also http://www.lx.it.pt/~mtf/GPSR
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob702.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'n','p','sigma'});
n           = getOption(parm,'n',       128);
p_spike     = getOption(parm,'p',      0.01);
sigma       = getOption(parm,'sigma',0.0005);
info.name   = 'blurspike';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',0); end;

% Set up the problem
data.signal          = sign(rand(n,n)-0.5).*(rand(n,n) <= p_spike);
data.noise           = sigma * randn(n*n,1);
data.op.Blur         = opBlur(n,n);
data.M               = data.op.Blur;
data.b               = data.M(reshape(data.signal,n*n,1),1) + data.noise;
data                 = completeOps(data);

% Additional information
info.title           = 'GPSR 2D Spikes';
info.thumb           = 'figProblem702';
info.citations       = {'FiguNowaWrig:2007'};
info.fig{1}.title    = 'Original binary spike image';
info.fig{1}.filename = 'figProblem702Image';
info.fig{2}.title    = 'Blurred image';
info.fig{2}.filename = 'figProblem702Blurred';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(data.signal), colormap gray;
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(reshape(data.b,n,n)), colormap gray;
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)
  
  if opts.update
     mn = min(data.b);
     mx = max(data.b);
     P = (data.b - mn) / (mx - mn);
     P = scaleImage(reshape(P,n,n),128,128);
     P = P(1:2:end,1:2:end,:);
     thumbwrite(P, info.thumb, opts);
  end
end
