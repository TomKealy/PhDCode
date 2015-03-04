function data = probSparseman(varargin)
%PROB052  Daubechies basis, "Man peeling orange" image.
%
%   PROB052 creates a problem structure.  The generated signal will
%   consist of N = 1024 by N grayscale image whose wavelet
%   coefficients are sparsified to retain only the P * N*N = 25,000
%   largest components. The image is randomly sampled at K = 96,000
%   points.
%
%   The following optional arguments are supported:
%
%   PROB052('k',K,'p',P,flags) is the same as above, but with a
%   K random measurements and preservation of a fraction P of the
%   largest wavelet coefficients. The 'noseed' flag can be specified
%   to suppress initialization of the random number generators. Each
%   of the parameter pairs and flags can be omitted.
%
%   Examples:
%   P = prob052;  % Creates the default 052 problem.
%
%   References:
%
%   [CandRomb:2006] E.J. Candes and J. Romberg, Sparsity and
%     incoherence in compressive sampling, Submitted for
%     publication, November, 2006.
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: probSparseman.m 1040 2008-06-26 20:29:02Z ewout78 $

[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'k','p'});
k           = getOption(parm,'k',96000); % Number of observations
p           = getOption(parm,'p',[]);    % Fraction of coefficients

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',3); end;

% Load the image and get wavelet coefficients
filename= sprintf('%sprob052_man.bmp', opts.datapath);
signal  = double(imread(filename)) / 256;
[m,n]   = size(signal);
wavelet = opWavelet(m,n,'Daubechies');

% Sparsify coefficients
p = 25000; if ~isempty(p), p = round(p*m*n); end
signal = sparsify(signal,p,wavelet);

% Set up the data
p = randperm(m*n); % Location of observations p(1:k)

% Set up the problem
data.x0              = wavelet(signal,2);
data.signal          = signal;
data.op.Daubechies   = wavelet;
data.op.Restrict     = opRestriction(m*n,p(1:k));
data.M               = data.op.Restrict;
data.B               = data.op.Daubechies;
data.b               = data.M(reshape(data.signal,m*n,1),1);
data                 = completeOps(data);

% Additional information
info.title           = 'Sparsified image of man with orange';
info.thumb           = 'figProblem052';
info.citations       = {'CandRomb:2006'};
info.fig{1}.title    = 'Complete image';
info.fig{1}.filename = 'figProblem052Full';
info.fig{2}.title    = 'Sampled points';
info.fig{2}.filename = 'figProblem052Samples';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  image(data.signal * 64), colormap gray;
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  P = max(0,min(1,data.signal)) * 64;
  P(p(k+1:end)) = 65;
  clrmap = [repmat(linspace(0,1,64)',1,3); 0.1 0.1 0.4];
  image(P), colormap(clrmap);
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)
  
  if opts.update
     P = scaleImage(data.signal,128,128);
     P = P(1:2:end,1:2:end,:);
     thumbwrite(P, info.thumb, opts);
  end
end
