function data = prob501(varargin)
%PROB501  Shepp-Logan phantom, partial Fourier along radial lines,
%         complex domain, total variation.
%
%   PROB501 creates a problem structure.  The generated signal will
%   consist of a N = 64 by N Shepp-Logan phantom. The signal is
%   measured along L = 22 radial lines in the frequency domain.
%
%   The following optional arguments are supported:
%
%   PROB501('l',L,'n',N,flags) is the same as above, but with a
%   phantom of size N by N, measured along L radial lines. The
%   'noseed' flag can be specified to suppress initialization of
%   the random number generators. Each of the parameter pairs and
%   flags can be omitted.
%
%   Examples:
%   P = prob501;  % Creates the default 501 problem.
%
%   References:
%
%   [CandRombTao:2006b] E.J. Candes, J. Romberg and T. Tao, Robust
%     uncertainty principles: exact signal reconstruction from
%     highly incomplete frequency information, IEEE Transactions on
%     Information Theory, 52:2 (2006), pp. 489-509.
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob501.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'l','n'});
n           = getOption(parm,'n',64); % Size of image
l           = getOption(parm,'l',22); % Number of radial lines
info.name   = 'phantom1';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',0); end;

% Set up the data
[M,Mh]  = RadialLines(n,l);
Mh(1,1) = 1;

% Set up the problem
data.signal          = ellipses(n);
data.mask            = Mh;
data.op.Restr        = opRestriction(n*n,find(data.mask));
data.op.FFT          = opFFT2d(n,n);
data.op.Haar         = opHaar2d(n,n);
data.M               = opFoG(data.op.Restr, data.op.FFT);
data.B               = data.op.Haar;
data.b               = data.M(reshape(data.signal,n*n,1),1);
data                 = completeOps(data);

% Additional information
info.title           = 'Shepp-Logan phantom';
info.thumb           = 'figProblem501';
info.citations       = {'CandRombTao:2006b'};
info.fig{1}.title    = 'The standard Shepp-Logan phantom';
info.fig{1}.filename = 'figProblem501SheppLogan';
info.fig{2}.title    = 'Radial lines in the Fourier domain';
info.fig{2}.filename = 'figProblem501RadialLinesFFT';
info.fig{3}.title    = 'Radial lines with FFT shift';
info.fig{3}.filename = 'figProblem501RadialLines';
info.fig{4}.title    = 'Half of the radial lines, shifted';
info.fig{4}.filename = 'figProblem501RadialLinesHalf';
info.fig{5}.title    = 'Sorted Haar wavelet coefficients';
info.fig{5}.filename = 'figProblem501SheppLoganHaarSort';
info.fig{6}.title    = 'Haar wavelet coefficients';
info.fig{6}.filename = 'figProblem501SheppLoganHaar';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  image(data.signal * 64), colormap gray;
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  spy(M,'b.');
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  spy(fftshift(M),'b.');
  updateFigure(opts, info.fig{3}.title, info.fig{3}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  spy(Mh,'b.');
  updateFigure(opts, info.fig{4}.title, info.fig{4}.filename)
  
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  H = data.op.Haar(reshape(data.signal,n*n,1),2);
  plot(sort(abs(H),'descend'),'b-'); xlim([1,length(H)]);
  updateFigure(opts, info.fig{5}.title, info.fig{5}.filename)
  
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(H,'b-'); xlim([1,length(H)]);
  updateFigure(opts, info.fig{6}.title, info.fig{6}.filename)
  
  if opts.update
     P = ellipses(256);
     P = (P(1:4:end,:) + P(2:4:end,:) + P(3:4:end,:) + P(4:4:end,:))/4;
     P = (P(:,1:4:end) + P(:,2:4:end) + P(:,3:4:end) + P(:,4:4:end))/4;
     thumbwrite(P, info.thumb, opts);
  end
end
