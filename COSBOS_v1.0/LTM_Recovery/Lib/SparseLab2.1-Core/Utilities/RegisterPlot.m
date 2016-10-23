function RegisterPlot(figname,nrflag)
% RegisterPlot -- Add legend with file name, date, flag
%  Usage
%    RegisterPlot(figname,nrflag)
%  Inputs
%    figname    string. In principle this can be any string.
%               In customary use, one supplies the name of the script
%               file generating the current display.
%    nrflag     string.  In principle this can be any string
%               In customary use, one supplies '{R} or '{NR}'
%               whether the plot is reproducible or not.
%
%  Side Effects
%    Near the bottom of the current figure, a legend is displayed
%    in a small font; the legend contains the figname and nrflag strings
%    and the current date.
%
%  See Also
%    AdaptDemo, BlockyDemo, IdealDemo, TourDemo
%
axes('Position',[.01 .01 .98 .05],'Visible','off');
h1 = text(.05,.5,figname); 
h2 = text(.80,.5,date);
h3 = text(.95,.5,nrflag);
set(h1,'FontSize',9);
set(h2,'FontSize',9);
set(h3,'FontSize',9);
    
    
    
%
% Copyright (c) 2006. David Donoho
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
