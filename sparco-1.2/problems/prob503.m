function data = prob503(varargin)
%PROB503  Shepp-Logan phantom, partial Fourier with sample mask,
%         complex domain, total variation.
%
%   PROB503 creates a problem structure.  The generated signal will
%   consist of a N = 256 by N Shepp-Logan phantom. The signal is
%   sampled at random locations in frequency domain generated
%   according to a probability density function.
%
%   The following optional arguments are supported:
%
%   PROB503('n',N,flags) is the same as above, but with a
%   phantom of size N by N. The 'noseed' flag can be specified to
%   suppress initialization of the random number generators. Both
%   the parameter pair and flags can be omitted.
%
%   Examples:
%   P = prob503;  % Creates the default 503 problem.
%
%   References:
%
%   [LustDonoPaul:2007] M. Lustig, D.L. Donoho and J.M. Pauly,
%     Sparse MRI: The application of compressed sensing for rapid MR
%     imaging, Submitted to Magnetic Resonance in Medicine, 2007.
%
%   [sparsemri] M. Lustig, SparseMRI,
%     http://www.stanford.edu/~mlustig/SparseMRI.html
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob503.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'n'});
n           = getOption(parm,'n',256);
info.name   = 'phantom2';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('twister',2000); end;

% Set up the data
pdf  = genPDF([n,n],5,0.33,2,0.1,0);
mask = genSampling(pdf,10,60);

% Set up the problem
data.signal          = ellipses(n);
data.noise           = randn(n)*0.01 + sqrt(-1)*randn(n)*0.01;
data.op.mask         = opMask(mask);
data.op.padding      = opPadding([n,n],[n,n]);
data.op.fft2d        = opFFT2C(n,n);
data.M               = opFoG(data.op.mask, data.op.padding, ...
                             data.op.fft2d);
data.b               = data.M(reshape(data.signal + data.noise,[n*n,1]),1);
data                 = completeOps(data);

% Additional information
info.title           = 'Shepp-Logan';
info.thumb           = 'figProblem503';
info.citations       = {'LustDonoPaul:2007','sparsemri'};
info.fig{1}.title    = 'Shepp-Logan phantom';
info.fig{1}.filename = 'figProblem503Phantom';
info.fig{2}.title    = 'Probability density function';
info.fig{2}.filename = 'figProblem503PDF';
info.fig{3}.title    = 'Sampling mask';
info.fig{3}.filename = 'figProblem503Mask';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(data.signal), colormap gray;
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(pdf), colormap gray;
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(mask), colormap gray
  updateFigure(opts, info.fig{3}.title, info.fig{3}.filename)
  
  if opts.update
     mn = min(min(data.signal + real(data.noise)));
     mx = max(max(data.signal + real(data.noise)));
     P = (data.signal + real(data.noise) - mn) / (mx - mn);
     P = scaleImage(P,128,128);
     P = P(1:2:end,1:2:end,:);
     thumbwrite(P, info.thumb, opts);
  end
end
