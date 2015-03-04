function data = prob603(varargin)
%PROB603  GPSR example: Smooth field.
%
%   PROB603 creates a problem structure.  The generated signal will
%   consist of an N = 64 by N image generated using the pwsmoothfield
%   function provided by earlier versions of GPSR. The signal is
%   measured using an M = 600 by N*N binary ensemble with unit norm
%   columns. Normally distributed noise with standard deviation
%   SIGMA = 0.001 is added to the measured signal.
%
%   The following optional arguments are supported:
%
%   PROB603('m',M,'n',N,'sigma',SIGMA,flags) is the same as above,
%   but with a measurement matrix of size M by N*N, a signal length
%   of N and a noise level of SIGMA. The 'noseed' flag can be
%   specified to suppress initialization of the random number
%   generators. Each of the parameter pairs and flags can be omitted.
%
%   Examples:
%   P = prob603;  % Creates the default 603 problem.
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
%   $Id: prob603.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'m','n','sigma'});
n           = getOption(parm,'n',       64); % Size of image (n x n)
m           = getOption(parm,'m',     1024); % Number of projections
sigma       = getOption(parm,'sigma',0.001); % Noise level
info.name   = 'yinyang';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',0); end;

% Set up the data
M = (2*round(rand(m,n*n))-1)/sqrt(n*n);

% Set up the problem
data.signal          = pwsmoothfield(n,1,0.05);
data.noise           = sigma * randn(m,1);
data.op.Matrix       = opMatrix(M);
data.op.Daubechies   = opWavelet(n,n,'Daubechies',2);
data.M               = data.op.Matrix;
data.B               = data.op.Daubechies;
data.b               = data.M(reshape(data.signal,n*n,1),1) + data.noise;
data                 = completeOps(data);

% Additional information
info.title           = 'GPSR 2D Compressed sensing';
info.thumb           = 'figProblem603';
info.citations       = {'FiguNowaWrig:2007'};
info.fig{1}.title    = 'Original image';
info.fig{1}.filename = 'figProblem603Image';
info.fig{2}.title    = 'Measurement matrix';
info.fig{2}.filename = 'figProblem603Matrix';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(data.signal), colormap gray;
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(M), colormap gray;
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)
  
  if opts.update
     P = scaleImage(data.signal,64,64);
     thumbwrite(P, info.thumb, opts);
  end
end
