function data = prob902(varargin)
%PROB902  Jittered sampling of cosine wave
%
%   PROB902 creates a problem structure.  The generated signal will
%   have length N = 1000 and consists of a cosine wave that is
%   sparse in the discrete cosine transform. The signal jitter
%   sampled; i.e. the sample locations are based on evenly spaced
%   intervals GAMMA apart, but are perturbed by random amounts in
%   the interval [-XI/2, XI/2], with default XI = GAMMA = 5.
%
%   The following optional argument calls are supported:
%
%   PROB902('k,',K,'n',N,'gamma',GAMMA,'xi',XI) is the same as
%   above but with a signal of length N having K non-zeros in the
%   discrete cosine domain. GAMMA gives the spacing of the regular
%   grid and XI prescribes the random perturbation as above.
%
%   Examples:
%   P = prob902;  % Creates the default 902 problem.
%
%   References:
%
%   [HennHerr:2007a] G. Hennenfent and F.J. Herrmann, Simply
%      denoise: wavefield reconstruction via coarse nonuniform
%      sampling, Technical Report, UBC Earth & Ocean Sciences,
%      August 2007.
%
%   [HennHerr:2007b] G. Hennenfent and F.J. Herrmann, Random
%      sampling: new insights into the reconstruction of
%      coarsely-sampled wavefields, SEG International Exposition
%      and 77th Annual Meeting, 2007.
%
%   [HennHerr:2007c] G. Hennenfent and F.J. Herrmann, Irregular
%      sampling: from aliasing to noise, EAGE 69th Conference &
%      Exhibition, 2007.
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Gilles Hennenfent, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob902.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'k','n','gamma','xi'});
k           = getOption(parm,'k',        3);
n           = getOption(parm,'n',     1000);
gamma       = getOption(parm,'gamma',    5);
xi          = getOption(parm,'xi',   gamma); % Jitter parameter
info.name   = 'jitter';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',1); rand('state',0); end;

% Set up the data
x0         = zeros(n,1);
p          = randperm(n);
x0(p(1:k)) = sign(randn(k,1)).*randn(k,1);
D          = opDCT(n);
signal     = D(x0,1);

% Set up the jitter sample indices
idx0    = 0:gamma:n-1; % Start at zero for modulo below
perturb = round(xi/2*rand(size(idx0)).*sign(randn(size(idx0))));
idx     = mod(idx0 + perturb,n)+1;

% Set up the problem
data.signal  = signal;
data.op.DCT  = D;
data.B       = data.op.DCT;
data.M       = opRestriction(n,idx);
data.b       = data.M(data.signal,1);
data.x0      = x0;
data         = completeOps(data);


% Additional information
info.title           = 'Jittered sampling of cosine wave';
info.thumb           = 'figProblem902';
info.citations       = {'HennHerr:2007a','HennHerr:2007b','HennHerr:2007c'};
info.fig{1}.title    = 'Signal';
info.fig{1}.filename = 'figProblem902Signal';
info.fig{2}.title    = 'Solution';
info.fig{2}.filename = 'figProblem902Solution';
info.fig{3}.title    = 'Back-projection';
info.fig{3}.filename = 'figProblem902BackProj';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(data.signal,'b-'); xlim([1,n]);
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(data.x0,'b-'); xlim([1,n]);
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(data.A(data.b,2),'b-');
  updateFigure(opts, info.fig{3}.title, info.fig{3}.filename)

 
  if opts.update
     P = ones(128,128,3);
     P = thumbPlot(P,1:n,data.A(data.b,2),[0,0,1]);
     P = (P(1:2:end,:,:) + P(2:2:end,:,:)) / 2;
     P = (P(:,1:2:end,:) + P(:,2:2:end,:)) / 2;
     thumbwrite(P, info.thumb, opts);
 end
end
