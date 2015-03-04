function data = prob703(varargin)
%PROB703  Curvelet basis, scratched fingerprint image.
%
%   PROB703 creates a problem structure.  The generated signal will
%   consist of an N = 128 by N grayscale image of a fingerprint with
%   data missing along straight lines.
%
%   The following optional arguments are supported:
%
%   PROB703('n',N,flags) is the same as above, but with an image size
%   of N by N. The 'noseed' flag can be specified to suppress
%   initialization of the random number generators. Both the
%   parameter pair and flags can be omitted.
%
%   Examples:
%   P = prob703;  % Creates the default 703 problem.
%
%   References:
%
%   [MaltMaioJainPrab:2003] D. Maltoni, D. Maio, A.K. Jain and
%     S. Prabhakar, Handbook of Fingerprint Recognition, Springer,
%     New York, 2003. See: http://bias.csr.unibo.it/fvc2000/
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob703.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'n'});
n           = getOption(parm,'n',128);
info.name   = 'finger';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',0); end;

% Set up the data
i = 4; j = 8;
x = double(imread(sprintf('%sprob703_FingerA1%02d_%d.tif',opts.datapath,i,j)));
x = scaleImage(x,n,n);

% Normalize image and apply color weighting
mn           = min(min(x));
mx           = max(max(x));
x            = 1- (x-mn) / (mx - mn);
x(x < 0.375) = 0;
x(x > 0.625) = 1;
idx          = find(x > 0 & x < 1);
x(idx)       = 0.5 + sin(4*pi*x(idx)) / 2;

% Set up the problem
data.signal      = x;
data.mask        = scratchMask(n,n,ceil(n/2));
data.op.Restrict = opRestriction(n*n,find(data.mask == 0));
data.op.Curvelet = opTranspose(opCurvelet2d(n,n));
data.M           = data.op.Restrict;
data.B           = data.op.Curvelet;
data.b           = data.M(reshape(data.signal,n*n,1),1);
data             = completeOps(data);

% Additional information
info.title           = 'Fingerprint on scratched surface';
info.thumb           = 'figProblem703';
info.citations       = {'MaltMaioJainPrab:2003'};
info.fig{1}.title    = 'Original fingerprint';
info.fig{1}.filename = 'figProblem703Fingerprint';
info.fig{2}.title    = 'Color reweighting';
info.fig{2}.filename = 'figProblem703Weights';
info.fig{3}.title    = 'Scratch pattern';
info.fig{3}.filename = 'figProblem703Mask';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  image((1-data.signal) * 64), colormap gray;
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  v = linspace(0,1,1024); w = sin(4*pi*v);
  w(v<(3/8.))   = -1;
  w(v > (5/8.)) =  1;
  w = w/2 + 0.5;
  plot(v,w,'b-'); xlim([0,1]);
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  spy(data.mask,'b.');
  updateFigure(opts, info.fig{3}.title, info.fig{3}.filename)
  
  if opts.update
     P = scaleImage(1 - data.signal,64,64);
     thumbwrite(P, info.thumb, opts);
  end
end
