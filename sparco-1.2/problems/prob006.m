function data = prob006(varargin)
%PROB006  Daubechies basis, Gaussian ensemble measurement basis,
%   piecewise cubic polynomial signal.
%
%   PROB006 creates a problem structure.  The generated signal will
%   consist of P = 5 cubic polynomial pieces and have length N = 2048.
%   The signal is measured using an M = 600 by N matrix with normally
%   distributed entries.
%
%   The following optional arguments are supported:
%
%   PROB006('m',M,'n',N,'p',P,flags) is the same as above, but with a
%   measurement matrix of size M by N and a signal length of N,
%   consisting of P cubic polynomial pieces. The 'noseed' flag can be
%   specified to suppress initialization of the random number
%   generators. Each of the parameter pairs and flags can be omitted.
%
%   Examples:
%   P = prob006;  % Creates the default 006 problem.
%
%   References:
%
%   [CandRomb:2004a] E.J. Candes and J. Romberg, Practical signal
%     recovery from random projections, Wavelet Applications in
%     Signal and Image Processing XI, Proc. SPIE Conf. 5914., 2004.
%
%   See also GENERATEPROBLEM.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: prob006.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse parameters and set problem name
[opts,varg] = parseDefaultOpts(varargin{:});
[parm,varg] = parseOptions(varg,{'noseed'},{'m','n','p'});
m           = getOption(parm,'m', 600); % Rows
n           = getOption(parm,'n',2048); % Columns
p           = getOption(parm,'p',   5); % Number of pieces
info.name   = 'p3poly';

% Return problem name if requested
if opts.getname, data = info.name; return; end;

% Initialize random number generators
if (~parm.noseed), randn('state',0); rand('state',0); end;

% Generate random piecewise cubic polynomial
idx = 1 + round(sort([0,rand(1,p-1),1]) * n);
signal = zeros(n,1);
for i=1:p
  x = linspace(-3,3,idx(i+1)-idx(i));
  c = round(randn(1,4) * 100) / 100;
  y = c(1) + c(2) * x + c(3) * x.^2 + c(4) * x.^3;
  signal(idx(i):idx(i+1)-1) = y;
  coefs{i} = c;
end

% Set up the problem
data.signal          = signal;
data.op.Daubechies   = opWavelet(n,1,'Daubechies');
data.op.Gaussian     = opGaussian(m,n,1);
data.M               = data.op.Gaussian;
data.B               = data.op.Daubechies;
data.b               = data.op.Gaussian(data.signal,1);
data.x0              = data.op.Daubechies(data.signal,2);
data                 = completeOps(data);

% Additional information
info.title           = 'Piecewise cubic polynomial';
info.thumb           = 'figProblem006';
info.citations       = {'CandRomb:2004a'};
info.fig{1}.title    = 'Piecewise polynomal function';
info.fig{1}.filename = 'figProblem006Signal';
info.fig{2}.title    = 'Wavelet coefficients of the piecewise polynomal function';
info.fig{2}.filename = 'figProblem006Coef';

% Plot figures
if opts.update || opts.show
  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(1:n,data.signal,'b-'); xlim([1,n]);
  updateFigure(opts, info.fig{1}.title, info.fig{1}.filename)

  figure(opts.figno); opts.figno = opts.figno + opts.figinc;
  plot(1:n,data.x0,'b-'); xlim([1,n]);
  updateFigure(opts, info.fig{2}.title, info.fig{2}.filename)
  
  if opts.update
     P = ones(128,128,3);
     P = thumbPlot(P,1:n,data.signal,[0,0,1]);
     P = (P(1:2:end,:,:) + P(2:2:end,:,:)) / 2;
     P = (P(:,1:2:end,:) + P(:,2:2:end,:)) / 2;
     thumbwrite(P, info.thumb, opts);
  end
end

% Add coefficient table
info.tab{1}.title    = ['Location and coefficients for the' ...
                        'cubic polynomial functions'];
info.tab{1}.filename = 'tablePiecePolyCoef';
info.tab{1}.celldata = {};
info.tab{1}.celltype = {};
info.tab{1}.celldata{1,1} = 'i';
info.tab{1}.celltype{1,1} = 'ml13';
info.tab{1}.celldata{2,1} = 'n(i)';
info.tab{1}.celltype{2,1} = 'ml';
for i=1:4
  if i==4, celltype = 'ml3'; else celltype = 'ml'; end;
    info.tab{1}.celldata{2+i,1} = sprintf('c_%d(i)',i-1);
    info.tab{1}.celltype{2+i,1} = celltype;
  end
  for i=1:p
    info.tab{1}.celldata{1,1+i} = sprintf('%d',i);
    info.tab{1}.celltype{1,1+i} = 'r13';
    info.tab{1}.celldata{2,1+i} = sprintf('%d',idx(i+1)-idx(i));
    info.tab{1}.celltype{2,1+i} = 'r';
  end
  for i=1:4
    if i==4, celltype = 'r3'; else celltype = 'r'; end;
    for j=1:p
      info.tab{1}.celldata{2+i,j+1} = sprintf('%5.3f', coefs{j}(i));
      info.tab{1}.celltype{2+i,j+1} = celltype;
    end
end

% Set the info field in data
data.info = info;
