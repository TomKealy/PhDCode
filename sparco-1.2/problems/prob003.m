function data = prob003(varargin)
%PROB003  DCT/Dirac dictionary, cosine/spikes signal.
%
%   PROB003 creates a problem structure.  The generated signal will
%   have length N = 1024 with K = 120 random spikes and C = 2 cosine
%   components.
%
%   The following optional arguments are supported:
%
%   PROB003('n',N,'k',K,'c',C,flags) is the same as above, but with a
%   signal length N, K random spikes and C cosine components. Set
%   C = [] to use the default cosine form
%
%   4.0 * cos(4*pi*t) + 2.0 * cos(12*pi*t)
%
%   for t = [0,1). The 'noseed' flag can be specified to suppress
%   initialization of the random number generators. Each of the
%   parameter pairs and flags can be omitted.
%
%   Examples:
%   P = prob003;  % Creates the default 003 problem.
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob003.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'c','k','n'});
c           = getOption(parm,'c',  []);
k           = getOption(parm,'k', 120);
n           = getOption(parm,'n',1024);
info.name   = 'cosspike';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',0); end;

% Cosine part of signal
t       = ((0:n-1) / n)';
odct    = opDCT(n);
cosCoef = zeros(n,1);
if isempty(c)
   cosCoef(5) = sqrt(n/2)*4;
   cosCoef(13)= sqrt(n/2)*2;
   % 4.0 * cos(4*pi*t) + 2.0 * cos(12*pi*t);
else
   p          = randperm(n); p = p(1:c);
   cosCoef(p) = randn(c,1) * sqrt(n/2);
end
cosine   = odct(cosCoef,1);

% Spike part of signal
p        = randperm(n); p = p(1:k);
spike    = zeros(n,1);
spike(p) = randn(k,1);


% Set up the problem
data.signal          = cosine + spike;
data.op.DCT          = opDCT(n);
data.op.Dirac        = opDirac(n);
data.op.Dict         = opDictionary(data.op.DCT, data.op.Dirac);
data.B               = data.op.Dict;
data.b               = data.signal;
data.x0              = [cosCoef; spike];
data                 = completeOps(data);

% Additional information
info.title           = 'Cosine with random spikes';
info.thumb           = 'figProblem003';
info.fig{1}.title    = 'Cosine with spikes';
info.fig{1}.filename = 'figProblem003SpikyCos';
info.fig{2}.title    = 'Cosine component of the signal';
info.fig{2}.filename = 'figProblem003Cosine';
info.fig{3}.title    = 'Spike component of the signal';
info.fig{3}.filename = 'figProblem003Spikes';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(t,data.signal,'b-'); xlim([0,1]);
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(t,cosine,'b-'); xlim([0,1]);
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(t,spike,'b-'); xlim([0,1]);
  updateFigure(opts, info.fig{3}.title, info.fig{3}.filename)
  
  if opts.update
     P = ones(128,128,3);
     P = thumbPlot(P,1:n,data.signal,[0,0,1]);
     P = (P(1:2:end,:,:) + P(2:2:end,:,:)) / 2;
     P = (P(:,1:2:end,:) + P(:,2:2:end,:)) / 2;
     thumbwrite(P, info.thumb, opts);
  end
end
