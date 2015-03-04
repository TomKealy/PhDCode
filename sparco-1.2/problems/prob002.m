function data = prob002(varargin)
%PROB002  Haar basis, block signal.
%
%   PROB002 creates a problem structure.  The generated signal will
%   have length N = 1024.
%
%   The following optional arguments are supported:
%
%   PROB002('n',N,flags) is the same as above, but with a signal length
%   N. The 'noseed' flag can be specified to suppress initialization of
%   the random number generators. Both the parameter pair and flags can
%   be omitted.
%
%   Example:
%   P = prob002;  % Creates the default 002 problem.
%
%   References:
%
%   [BuckDono:1995] J. Buckheit and D.L. Donoho, WaveLab and
%     reproducible research, in: A. Antoniadis (Ed.), Wavelets and
%     Statistics, Springer, Berlin, 1995. 
%
%   [ChenDonoSaun:1998] S.S. Chen and D.L. Donoho and M.A. Saunders,
%     Atomic Decomposition by Basis Pursuit, SIAM Journal on
%     Scientific Computing 20:1 (1998), pp. 33-61.
%
%   [DonoJohn:1993] D.L. Donoho and I.M. Johnstone, Ideal spatial
%     adaptation by wavelet shrinkage, Biometrika, 81:3 (1994),
%     pp. 425-455.
%
% See also GENERATEPROBLEM.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob002.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'n'});
n           = getOption(parm,'n',1024);
info.name   = 'blocksig';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Set up the problem
data.signal          = makesig('Blocks',n)';
data.op.Haar         = opHaar(n);
data.B               = data.op.Haar;
data.b               = data.signal;
data.x0              = data.op.Haar(data.signal,2);
data                 = completeOps(data);

% Additional information
info.title           = 'Block signal';
info.thumb           = 'figProblem002';
info.citations       = {'BuckDono:1995','ChenDonoSaun:1998','DonoJohn:1993'};
info.fig{1}.title    = 'Blocks signal';
info.fig{1}.filename = 'figProblem002Blocks';
info.fig{2}.title    = 'Haar coefficients of the blocks signal';
info.fig{2}.filename = 'figProblem002Coef';

% Set the info field in data
data.info = info;

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot((1:n)./n, data.signal,'b-');
  xlim([0,1]); ylim([-3,6]);
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)
       
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(data.op.Haar(data.signal,2),'b-');
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)

  if opts.update
     P = ones(128,128,3);
     P = thumbPlot(P,1:n,data.signal,[0,0,1]);
     P = (P(1:2:end,:,:) + P(2:2:end,:,:)) / 2;
     P = (P(:,1:2:end,:) + P(:,2:2:end,:)) / 2;
     thumbwrite(P, info.thumb, opts);
  end
end
