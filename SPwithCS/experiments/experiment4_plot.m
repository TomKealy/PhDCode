%% experiment4_plot.m - Plots Figure 7
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

load('experiment4_data5');    
signalNorm5 = signalNorm;
errorNorm5 = errorNorm;

load('experiment4_data10');    
signalNorm10 = signalNorm;
errorNorm10 = errorNorm;

load('experiment4_data15');    
signalNorm15 = signalNorm;
errorNorm15 = errorNorm;

N = 1024*4;
B = 1/256;
mmmax = 13;
kkmax = 50;

snr5 = zeros(kkmax,mmmax);
snr10 = zeros(kkmax,mmmax);
snr15 = zeros(kkmax,mmmax);

for kk=1:kkmax,
    for mm = 1:mmmax,
        snr5(kk,mm) = 20*log10(signalNorm5(kk)/errorNorm5(mm,kk));
        snr10(kk,mm) = 20*log10(signalNorm10(kk)/errorNorm10(mm,kk));
        snr15(kk,mm) = 20*log10(signalNorm15(kk)/errorNorm15(mm,kk));
    end
end

M_vec5 = zeros(mmmax,1);
M_vec10 = zeros(mmmax,1);
M_vec15 = zeros(mmmax,1);
for mm=1:mmmax,
    M_vec5(mm) = N*B*5*(mm+1)/2;
    M_vec10(mm) = N*B*10*(mm+1)/2;
    M_vec15(mm) = N*B*15*(mm+1)/2;
end

snr1_quant = quantile(snr5,[.1]);
snr2_quant = quantile(snr10,[.1]);
snr3_quant = quantile(snr15,[.1]);

figure(1)

plot(M_vec5/(N*B*5),snr1_quant,'-o','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
hold on
plot(M_vec10/(N*B*10),snr2_quant,'r:s','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
plot(M_vec15/(N*B*15),snr3_quant,'g--^','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])

axis([1 M_vec15(end)/(N*B*15) 0 250])

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
hLegend = legend('$K = 5$','$K = 10$','$K = 15$','Location','SouthEast');
set(hLegend,'LineWidth',2,'Color','w','Interpreter','LaTex');

set(gca, 'Position', positiondata);
export_fig 'experiment4a.pdf' -nocrop

figure(2)

plot(M_vec5,snr1_quant,'-o','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
hold on
plot(M_vec10,snr2_quant,'r:s','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
plot(M_vec15,snr3_quant,'g--^','LineWidth',2,'MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])

axis([1 M_vec15(end) 0 250])

xlabel('Measurements ($M$)','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
ylabel('Recovery SNR (dB)','FontSize',fontsize,'FontName','Times New Roman','Interpreter','LaTex');
set(gca, ...
    'LineWidth' , 2                     , ...
    'FontSize'  , fontsize              , ...
    'FontName'  , 'Times New Roman'     , ...
    'YTick'     , [0 50 100 150 200 250]      , ...
    'YGrid'     , 'on'                  , ...
    'Color'     , 'none'                , ...
    'Units'     , 'inches'              );
set(gcf, 'Color', 'none');
hLegend = legend('$K = 5$','$K = 10$','$K = 15$','Location','SouthEast');
set(hLegend,'LineWidth',2,'Color','w','Interpreter','LaTex');

set(gca, 'Position', positiondata);
export_fig 'experiment4b.pdf' -nocrop
