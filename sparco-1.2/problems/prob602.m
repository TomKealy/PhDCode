function data = prob602(varargin)
%PROB602  2D Haar basis, binary measurement basis, b/w image.
%
%   PROB602 creates a problem structure. The generated signal will
%   consist of a 64 by 64 black and white image of a soccer ball.
%   The signal is measured using an M by 4096 binary ensemble, with
%   M = 3200.
%
%   The following optional arguments are supported:
%
%   PROB602('m',M,'scale',SCALE,flags) is the same as above,
%   but with a an image size of 32 by 32 if SCALE = 1 and 64 by 64 if
%   SCALE = 2. The measurement matrix has a size M by (32*SCALE)^2.
%   Leaving M unspecified gives a measurement matrix of size
%   800*SCALE^2 by (32*SCALE)^2. The 'noseed' flag can be specified
%   to suppress initialization of the random number generators. Each
%   of the parameter pairs and flags can be omitted.
%
%   Examples:
%   P = prob602;  % Creates the default 602 problem.
%
%   References:
%
%   [TakhLaskWakiEtAl:2006] D. Takhar and J. N. Laska and M. Wakin
%     and M. Duarte and D. Baron and S. Sarvotham and K. K. Kelly
%     and R. G. Baraniuk, A new camera architecture based on
%     optical-domain compression, in Proceedings of the IS&T/SPIE
%     Symposium on Electronic Imaging: Computational Imaging,
%     January 2006, Vol. 6065.
%
%   The soccerball image is taken from Wikipedia.
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

% scale  1 is the 32 pixel image, 2 (default) is the 64 pixel image.


%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob602.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'m','scale'});
scale       = getOption(parm,'scale', 2); % 1 or 2
m           = getOption(parm,'m',    []);
info.name   = 'soccer2';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',0); end;

% Set up the data
n = 32 * scale;
if isempty(m), m = 800 * scale^2; end;

% Set up the problem
data.signal      = double(imread(sprintf('%sprob602_ball%d.png',opts.datapath,n)));
data.op.Binary   = opBinary(m,n*n);
data.op.Haar     = opHaar2d(n,n);
data.M           = data.op.Binary;
data.B           = data.op.Haar;
data.b           = data.M(reshape(data.signal,n*n,1),1);
data.x0          = [];
data             = completeOps(data);

% Additional information
info.title           = 'Soccer ball';
info.thumb           = 'figProblem602';
info.citations       = {'TakhLaskWakiEtAl:2006'};
info.fig{1}.title    = 'Original signal';
info.fig{1}.filename = 'figProblem602Signal';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  imagesc(data.signal); colormap gray;
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)
  
  if opts.update
     P = double(imread([opts.datapath,'prob602_ball64.png']));
     thumbwrite(P, info.thumb, opts);
  end
end
