function data = prob701(varargin)
%PROB701  GPSR example: Daubechies basis, blurred Photographer.
%
%   PROB701 creates a problem structure.  The generated signal will
%   consist of the 256 by 256 grayscale 'photographer' image. The
%   signal is blurred by convolution with an 8 by 8 blurring mask and
%   normally distributed noise with standard deviation SIGMA = 0.0055
%   is added to the final signal.
%
%   The following optional arguments are supported:
%
%   PROB701('sigma',SIGMA,flags) is the same as above, but with the
%   noise level set to SIGMA. The 'noseed' flag can be specified to
%   suppress initialization of the random number generators. Both the
%   parameter pair and flags can be omitted.
%
%   Examples:
%   P = prob701;  % Creates the default 701 problem.
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
%   $Id: prob701.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'sigma'});
sigma       = getOption(parm,'sigma', sqrt(2) / 256);
info.name   = 'blurrycam';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); end;

% Set up the data
signal = imread(sprintf('%sprob701_Camera.tif', opts.datapath));
[m,n]  = size(signal);

% Set up the problem
data.signal          = double(signal) / 256;
data.op.Blur         = opBlur(m,n);
data.op.Daubechies   = opWavelet(m,n,'Daubechies',2);
data.M               = data.op.Blur;
data.B               = data.op.Daubechies;
data.b               = data.M(reshape(data.signal,m*n,1),1);
data.b               = data.b + sigma * randn(m*n,1);
data                 = completeOps(data);

% Additional information
info.title           = 'GPSR 2D Deblur';
info.thumb           = 'figProblem701';
info.citations       = {'FiguNowaWrig:2007'};
info.fig{1}.title    = 'Original image';
info.fig{1}.filename = 'figProblem701Image';
info.fig{2}.title    = 'Blurred image';
info.fig{2}.filename = 'figProblem701Blurred';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(data.signal), colormap gray;
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(reshape(data.b,m,n)), colormap gray;
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)
  
  if opts.update
     P =  scaleImage(reshape(data.b,m,n),64,64);
     thumbwrite(P, info.thumb, opts);
  end
end
