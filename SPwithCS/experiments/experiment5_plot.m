%% experiment5_plot.m - Plots Figure 8
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

load('experiment5_data1');    

N = 1024*4;
B = 1/256;
K = 5;

mmmax = 13;
M_max = N*B*K*mmmax;

dpssSizeMin = N*B; 
dpssSizeMax = 38;

kkmax = 50;

snr1 = zeros(kkmax,mmmax);
snr2 = zeros(kkmax,mmmax);
snr3 = zeros(kkmax,mmmax);

for kk=1:kkmax,
    for mm = 1:mmmax,
        snr1(kk,mm) = 20*log10(signalNorm(kk)/errorNorm1(mm,kk));
        snr2(kk,mm) = 20*log10(signalNorm(kk)/errorNorm2(mm,kk));
        snr3(kk,mm) = 20*log10(signalNorm(kk)/errorNorm3(mm,kk));
    end
end

M_vec = zeros(mmmax,1);
for mm=1:mmmax,
    M_vec(mm) = N*B*K*(mm+1)/2;
end

figure(1)

snr1_quant = quantile(snr1,[.1]);
snr2_quant = quantile(snr2,[.1]);
snr3_quant = quantile(snr3,[.1]);

plot(M_vec/(N*B*K),snr1_quant,'-o','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
hold on
plot(M_vec/(N*B*K),snr2_quant,'r:s','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
plot(M_vec/(N*B*K),snr3_quant,'g--^','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])

axis([1 M_vec(end)/(N*B*K) 0 250])

set(gca,'Units','inches');
positiondata = get(gca,'Position');
positiondata(1) = 1;
positiondata(2) = .64;

xlabel('Oversampling factor $\left(\frac{M}{2NWK}\right)$','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Recovery SNR (dB)','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
set(gca, ...
    'LineWidth' , 2                     , ...
    'FontSize'  , fontsize              , ...
    'FontName'  , 'Times New Roman'     , ...
    'YTick'     , [0 50 100 150 200 250]      , ...
    'YGrid'     , 'on'                  , ...
    'Color'     , 'none'                );
set(gcf, 'Color', 'none');
hLegend = legend('Gaussian matrix','Random demodulator','Random sampling','Location','SouthEast');
set(hLegend,'LineWidth',2,'Color','w','Interpreter','LaTex','FontSize',fontsize2);

set(gca, 'Position', positiondata);
export_fig 'experiment5.pdf' -nocrop
