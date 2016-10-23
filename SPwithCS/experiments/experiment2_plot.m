%% experiment2_plot.m - Plots Figure 6
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

load experiment2_data1;

N = 1024*4;
B = 1/256;
dpssSizeInit = N*B-8;
jjmax = 14;

snr1 = 20*log10(repmat(signalNorm,jjmax,1)./errorNorm1);
snr2 = 20*log10(repmat(signalNorm,jjmax,1)./errorNorm2);
snr3 = 20*log10(repmat(signalNorm,jjmax,1)./errorNorm3);

snr1_quant = quantile(snr1',.1);
snr2_quant = quantile(snr2',.1);
snr3_quant = quantile(snr3',.1);

figure(1)

plot([dpssSizeInit-2:2:dpssSizeInit+2*jjmax-4],snr3_quant,'-o','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
hold on
plot([dpssSizeInit-2:2:dpssSizeInit+2*jjmax-4],snr2_quant,'r:s','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
plot([dpssSizeInit-2:2:dpssSizeInit+2*jjmax-4],snr1_quant,'g--^','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])

axis([6 32 0 115])

set(gca,'Units','inches');
positiondata = get(gca,'Position');
positiondata(1) = 1;
positiondata(2) = .64;

xlabel('DPSS vectors per band ($k$)','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Recovery SNR (dB)','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
set(gca, ...
    'LineWidth' , 2                     , ...
    'FontSize'  , fontsize              , ...
    'FontName'  , 'Times New Roman'     , ...
    'YTick'     , [0 25 50 75 100]      , ...
    'YGrid'     , 'on'                  , ...
    'XTick'     , [8 12 16 20 24 28 32] , ...
    'Color'     , 'none'                );
set(gcf, 'Color', 'none');
hLegend = legend('MSNR = 100dB','MSNR = 50dB','MSNR = 25dB','Location','NorthWest');
set(hLegend,'LineWidth',2,'Color','w','Interpreter','LaTex');

set(gca, 'Position', positiondata);
export_fig 'experiment2.pdf' -nocrop
