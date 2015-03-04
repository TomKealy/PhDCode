function data = prob502(varargin)
%PROB502  Angiogram, partial Fourier with sample mask, complex domain,
%   total variation.
%
%   PROB502 creates a problem structure.  The generated signal will
%   consist of an N = 100 by N synthetic angiogram image. The signal
%   is sampled at random locations in frequency domain generated
%   according to a probability density function.
%
%   The following optional arguments are supported:
%
%   PROB502('n',N,flags) is the same as above, but with an angiogram
%   of size N by N. The 'noseed' flag can be specified to suppress
%   initialization of the random number generators. Both of the
%   parameter pair and flags can be omitted.
%
%   Examples:
%   P = prob502;  % Creates the default 502 problem.
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
%   $Id: prob502.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'n'});
n           = getOption(parm,'n',100);
info.name   = 'angiogram';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',0); end;

% Set up the data
if (n == 100)
   % Load default angiogram, pdf and mask; this avoid having
   % to use the 'twister' method for rand and makes the code
   % compatible with older versions of Matlab.
   x = load(sprintf('%sprob502_angiogram.mat', opts.datapath));
   signal = x.angiogram;
   pdf    = x.pdf;
   mask   = x.mask;
   
   % ----------------------------------------------------------
   % Note: the above data is generated using the following code
   % ----------------------------------------------------------
   % randn('state',0); rand('twister',16000);
   % pdf    = genPDF([n,n],5,0.33,2,0.1,0);
   % mask   = genSampling(pdf,10,60);
   % signal = generateAngiogram(n,n);
   % ----------------------------------------------------------
else
   % Generate the data
   pdf    = genPDF([n,n],5,0.33,2,0.1,0);
   mask   = genSampling(pdf,10,60);
   signal = generateAngiogram(n,n);
end

% Set up the problem
data.pdf             = pdf;
data.mask            = mask;
data.signal          = signal;
data.op.mask         = opMask(mask);
data.op.padding      = opPadding([n,n],[n,n]);
data.op.fft2d        = opFFT2C(n,n);
data.M               = opFoG(data.op.mask, data.op.padding, ...
                             data.op.fft2d);
data.b               = data.M(reshape(data.signal,[n*n,1]),1);
data                 = completeOps(data);

% Additional information
info.title           = 'Artificial angiogram';
info.thumb           = 'figProblem502';
info.citations       = {'LustDonoPaul:2007','sparsemri'};
info.fig{1}.title    = 'Artificial angiogram';
info.fig{1}.filename = 'figProblem502Image';
info.fig{2}.title    = 'Probability density function';
info.fig{2}.filename = 'figProblem502PDF';
info.fig{3}.title    = 'Sampling mask';
info.fig{3}.filename = 'figProblem502Mask';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(data.signal), colormap gray;
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)
  
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(data.pdf), colormap gray;
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(data.mask), colormap gray;
  updateFigure(opts, info.fig{3}.title, info.fig{3}.filename)

  if opts.update
     mn = min(min(data.signal));
     mx = max(max(data.signal));
     P = (data.signal - mn) / (mx - mn);
     P = scaleImage(P,128,128);
     P = P(1:2:end,1:2:end,:);
     thumbwrite(P, info.thumb, opts);
  end
end
