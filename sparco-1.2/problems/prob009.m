function data = prob009(varargin)
%PROB009  Heaviside basis, block signal.
%
%   PROB009 creates a problem structure.  The generated signal will
%   have length N = 128.
%
%   The following optional argument calls are supported:
%
%   PROB009('n',N) is the same as above but with a signal length N.
%
%   Examples:
%   P = prob009;  % Creates the default 009 problem.
%
%   References:
%
%   [ChenDonoSaun:1998] S.S. Chen and D.L. Donoho and M.A. Saunders,
%     Atomic Decomposition by Basis Pursuit, SIAM Journal on
%     Scientific Computing 20:1 (1998), pp. 33-61.
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob009.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{},{'n'});
n           = getOption(parm,'n',128);
info.name   = 'blkheavi';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Set up the problem
data.signal          = makesig('Blocks',n)';
data.op.Heaviside    = opHeaviside(n);
data.B               = data.op.Heaviside;
data.b               = data.signal;
data.x0              = [data.signal(1); data.signal(2:end)-data.signal(1:end-1)];
data                 = completeOps(data);

% Additional information
info.title           = 'Block signal with Heaviside matrix';
info.thumb           = 'figProblem009';
info.citations       = {'ChenDonoSaun:1998','DonoJohn:1993'};
info.fig{1}.title    = 'Signal';
info.fig{1}.filename = 'figProblem009Signal';
info.fig{2}.title    = 'Solution';
info.fig{2}.filename = 'figProblem009Solution';

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
  
  if opts.update
     P = thumbFromOp(data.A,64,64,16,16,1); % Greyscale
     thumbwrite(P, info.thumb, opts);
  end
end
