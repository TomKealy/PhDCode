%% experiment1_plot.m - Plots Figure 5
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

clear all;
close all;
addpath('export_fig')

fontsize = 18;
fontsize2 = 16;

load experiment1_data1

N = 1024*4;
B = 1/256;
dpssSizeInit = N*B-8;
jjmax = 20;

snr1 = 20*log10(repmat(signalNorm,jjmax,1)./errorNorm1);
snr2 = 20*log10(repmat(signalNorm,jjmax,1)./errorNorm2);

figure(1)
snr1_quant = quantile(snr1',.1);
snr2_quant = quantile(snr2',.1);

plot([dpssSizeInit-2:2:dpssSizeInit+2*jjmax-4],snr1_quant,'-o','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
hold on
plot([dpssSizeInit-2:2:dpssSizeInit+2*jjmax-4],snr2_quant,'r:s','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])

axis([6 44 0 250])
set(gca,'Units','inches');
positiondata = get(gca,'Position');
positiondata(1) = 1;
positiondata(2) = .625;

xlabel('DPSS vectors per band ($k$)','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Recovery SNR (dB)','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
set(gca, ...
    'LineWidth'  , 2                      , ...
    'FontSize'   , fontsize               , ...
    'FontName'   , 'Times New Roman'      , ...
    'YTick'      , [0 50 100 150 200 250] , ...
    'YGrid'      , 'on'                   , ...
    'XTick'      , [8 16 24 32 40]        , ...
    'Color'      , 'none'                 );
set(gcf, 'Color', 'none');
hLegend = legend('BBCoSaMP($A,\Psi,y,K$)','BBCoSaMP($A \Psi,I,y,K$)','Location','NorthWest');
set(hLegend,'LineWidth',2,'Interpreter','LaTex','FontSize',fontsize2);

set(gca, 'Position', positiondata);
export_fig 'experiment1a.pdf' -nocrop

figure(2)
pe = (sign(snr1 - snr2 - .1*snr1)+1)/2;
plot([dpssSizeInit-2:2:dpssSizeInit+2*jjmax-4],1-mean(pe,2),'-o','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])

axis([6 44 0 1])

xlabel('DPSS vectors per band ($k$)','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Probability','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
set(gca, ...
    'LineWidth'  , 2                      , ...
    'FontSize'   , fontsize               , ...
    'FontName'   , 'Times New Roman'      , ...
    'YTick'      , [0 .25 .5 .75 1]       , ...
    'XTick'      , [8 16 24 32 40]        , ...
    'Color'      , 'none'                 , ...
    'Units'      , 'inches'               );
set(gcf, 'Color', 'none');

set(gca, 'Position', positiondata);
export_fig 'experiment1b.pdf' -nocrop


figure(3)

[~,V] = dpss(N,N*B/2,44);

V_trim = V(6:44);
semilogy([6:44],max(V_trim,eps),'-o','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0]);

axis([6 44 10^(-16)/2 1])

xlabel('DPSS vectors per band ($k$)','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('$\lambda_{N,W}^{(k+1)}$','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex','rot',0);
set(gca, ...
    'LineWidth'  , 2                      , ...
    'FontSize'   , fontsize               , ...
    'FontName'   , 'Times New Roman'      , ...
    'YMinorTick' , 'on'                   , ...
    'YTick'      , [10^(-15) 10^(-10) 10^(-5) 10^(0)] , ...
    'XTick'      , [8 16 24 32 40]        , ...
    'Color'      , 'none'                 , ...
    'Units'      , 'inches'               );
set(gcf, 'Color', 'none');

set(gca, 'Position', positiondata);
export_fig 'experiment1c.pdf' -nocrop
