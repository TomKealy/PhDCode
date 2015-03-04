function data = prob010(varargin)
%PROB010  Normalized Heaviside (unit columns), block signal.
%
%   PROB010 creates a problem structure.  The generated signal will
%   have length N = 1024.
%
%   The following optional argument calls are supported:
%
%   PROB010('n',N) is the same as above but with a signal length N.
%
%   Examples:
%   P = prob010;  % Creates the default 010 problem.
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
%   $Id: prob010.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{},{'n'});
n           = getOption(parm,'n',1024);
info.name   = 'blknheavi';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Set up the data
signal = makesig('Blocks',n)';

% Set up the problem
data.signal       = signal;
data.op.Heaviside = opHeaviside(n,1); % Scale columns
data.B            = data.op.Heaviside;
data.b            = data.signal;
data.x0           = [signal(1); signal(2:end)-signal(1:end-1)];
data.x0           = data.x0 .* sqrt(n:-1:1)';
data              = completeOps(data);

% Additional information
info.title           = 'Block signal with normalized Heaviside matrix';
info.thumb           = 'figProblem010';
info.citations       = {'ChenDonoSaun:1998','DonoJohn:1993'};
info.fig{1}.title    = 'Signal';
info.fig{1}.filename = 'figProblem010Signal';
info.fig{2}.title    = 'Solution';
info.fig{2}.filename = 'figProblem010Solution';

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
     op = opHeaviside(18,1);
     P = thumbFromOp(op,64,64,16,16,1); % Greyscale
     thumbwrite(P, info.thumb, opts);
  end
end
