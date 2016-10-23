%% experiment7_plot.m - Plots Figure 9b
%
% Most recent change - 9/15/2011
%
% Copyright 2011, Mark Davenport, Michael Wakin
%
% This file is part of DPSS Approximation and Recovery Toolbox version 1.0.
%
%    DART is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    DART is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DART.  If not, see <http://www.gnu.org/licenses/>.

close all; clear all;
addpath('../');
addpath('export_fig')

fontsize = 18;

s = RandStream.create('mt19937ar','seed',5489);
RandStream.setDefaultStream(s);

N = 1024*4;
B = 1/256;
K = 5;
M = 4*N*B*K;
S = 17*K;

dpssSize = round(5.5*M/(N*B*K)+5);
dpssWidth = N*B/2;

senseOpts.sampler = 'rdmdemod';
senseOpts.pnseq = sign(randn(N,1)); 
senseOpts.subsamp = N/M;

sigOpts.N = N; sigOpts.B = B; sigOpts.K = K;
sigOpts.type = 'randtones';  
sigOpts.TPB = 50; 

cosampOpts.N = N;
cosampOpts.maxiter = 100;
cosampOpts.tol = 1e-3;
cosampOpts.PsiRecover = dpss(N,dpssWidth,dpssSize); 
cosampOpts.modulator = generateModulator(N,B); 
cosampOpts.verbose = 0;
cosampOpts.detectAlg = 'blockOMP';
cosampOpts.debias = 1;

[x supp] = generateSignal(sigOpts); 
x = x/norm(x);

cosampOpts.normBound = 1+1e-13;

y = sense(x,0,senseOpts); 

[x_dpss, ~, ~] = dpssCosamp(y,senseOpts,K,B,cosampOpts);
[x_fft,~] = fftCosamp(y,senseOpts,S,cosampOpts);

X = abs(fftshift(fft(x)))/sqrt(N);
X_dpss = abs(fftshift(fft(x_dpss)))/sqrt(N);
X_fft = max(abs(fftshift(fft(x_fft)))/sqrt(N),10e-6);

figure(1)
semilogy(linspace(-.5,.5,N),X,'b','Linewidth',2);
hold on
semilogy(linspace(-.5,.5,N),X_dpss,'r:','Linewidth',2)
semilogy(linspace(-.5,.5,N),X_fft,'go','Linestyle','none','MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 1 0])

set(gca,'Units','inches');
positiondata = get(gca,'Position');
positiondata(1) = 1;
positiondata(2) = .64;
positiondata(3) = .99*positiondata(3);
positiondata(4) = .99*positiondata(4);

set(gca, ...
    'LineWidth'  , 2                      , ...
    'FontSize'   , fontsize               , ...
    'FontName'   , 'Times New Roman'      , ...
    'XLim'       , [-.5 .5]               , ...
    'XTick'      , [-.5 -.25 0 .25 .5]    , ...
    'YLim'       , [10e-6 1]              , ...
    'YTick'      , [1e-5 1e-4 1e-3 1e-2 1e-1 1]  , ...
    'Color'      , 'none'                  );
set(gcf, 'Color', 'none');

hYlabel = ylabel('$|\widehat{X}(f)|$','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex','rot',0,'HorizontalAlignment','right');
hXlabel = xlabel('$f$','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex'); 

xl = get(gca,'xlabel');
yl = get(gca,'ylabel');
orig_label_units = get( xl, 'Units');
set([xl yl], 'Units','pixels');
xp = get( xl, 'Position');
yp = get( yl, 'Position'); 

[hx,hy] = format_ticks(gca,{'$-\frac{1}{2}$','$-\frac{1}{4}$','$0$','$\frac{1}{4}$','$\frac{1}{2}$'}, ...
                           [] ,...
                           [-.5 -.25 0 .25 .5], ...
                           [], ...
                           [],[],[],'FontSize',fontsize,'FontName','Times New Roman');
                       
xp(1) = 260; 
yp(1) = -30; 
yp(2) = 195; 

set( xl, 'Position', xp);
set( yl, 'Position', yp);
set([xl yl], 'Units', orig_label_units);

hLegend = legend('Original','DPSS','DFT','Location','SouthEast');
set(hLegend,'LineWidth',2,'Color','w','Interpreter','LaTex');

set(gca, 'Position', positiondata);
export_fig 'experiment7.pdf' -nocrop

disp(['SNR1: ' num2str(20*log10(norm(x)/norm(x-x_dpss))) ' -- SNR2: ' num2str(20*log10(norm(x)/norm(x-x_fft)))])