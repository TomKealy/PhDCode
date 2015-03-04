function updateFigure(opts, figTitle, figFilename)

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: updateFigure.m 1040 2008-06-26 20:29:02Z ewout78 $

% Ensure default values are available
opts.linewidth  = getOption(opts,'linewidth', []);
opts.fontsize   = getOption(opts,'fontsize',  []);
opts.markersize = getOption(opts,'markersize',[]);

% Output the plots
if opts.update
  % Set the line width, font size and marker size
  chld = [gca; get(gca,'Children')];
  lnwd = ones(length(chld),1) * NaN;
  fnts = ones(length(chld),1) * NaN;
  mrks = ones(length(chld),1) * NaN;
  for i=1:length(chld)
    conf = get(chld(i));
    if ~isempty(opts.linewidth) && isfield(conf,'LineWidth')
      lnwd(i) = get(chld(i),'LineWidth');
      if (lnwd(i) == 0.5) % Default
        set(chld(i),'Linewidth',opts.linewidth);
      end
    end
    if ~isempty(opts.fontsize) && isfield(conf,'FontSize')
      fnts(i) = get(chld(i),'FontSize');
      if (fnts(i) == 10) % Default
        set(chld(i),'FontSize',opts.fontsize);
      end
    end
    if ~isempty(opts.markersize) && isfield(conf,'MarkerSize')
      mrks(i) = get(chld(i),'MarkerSize');
      if (mrks(i) == 6) % Default
        set(chld(i),'MarkerSize',opts.markersize);
      end
    end
  end
  
  for i=1:length(opts.figtype)
     updateFigureType(opts.update, 0, opts.figtype{i}, ...
                      opts.figpath, figTitle, figFilename);
  end

  % Restore the line-widths, font size
  for i=1:length(chld)
    if ~isnan(lnwd(i))
      set(chld(i),'LineWidth',lnwd(i));
    end
    if ~isnan(fnts(i))
      set(chld(i),'FontSize',fnts(i));
    end
    if ~isnan(mrks(i))
      set(chld(i),'MarkerSize',mrks(i));
    end
  end
    
end

% Show the plot
if opts.show
  updateFigureType(0,opts.show,'','',figTitle,'');
end



function updateFigureType(update,show,figtype,figpath,figTitle,figFilename)
filename = [figpath,figFilename];

switch lower(figtype)
 case {'pdf'}
   cmdPostprocess = sprintf('!pdfcrop %s.pdf %s.pdf >& /dev/null', ...
                            filename, filename);
 otherwise
   cmdPostprocess = [];
end

[figtype,figext] = getFigureExt(figtype);

% Print the figure for output (without title)
if update
  evalc(sprintf('print -d%s %s.%s;', figtype, filename, figext));
  
  if ~isempty(cmdPostprocess)
    eval(cmdPostprocess);
  end
end

% Add title if needed
if show
  title(figTitle);
end
