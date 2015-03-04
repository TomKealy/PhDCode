% Figures for paper
%   From Theory to Practice: 
%     Sub-Nyquist Sampling of Sparse Wideband Analog  Signals
%  By Moshe Mishali and Yonina C. Eldar

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main figure - support recovery for #channels
load Results_100Shift1
Data = squeeze(mean(Results(1:SigInd,:,:),1));
FigNum=1;
h=figure(FigNum);
axes('Parent',h,'Position',[0.1369 0.1475 0.7681 0.7775],...
    'LineWidth',1.5,...
    'FontSize',9);
image(Channel_axis,SNR_axis,256*Data')
colormap(gray(256))
colorbar('YTickLabel',{'0.2','0.4','0.6','0.8','1'})
xlabel('\# sampling channels ($\bar{m}$)','Interpreter','latex','FontSize',10);
ylabel('SNR (dB)','Interpreter','latex','FontSize',10,'Margin',2);
ax=get(h,'CurrentAxes');
set(ax,'YDir','normal')
axis([20 100 -20 30])




%%%%%%%%%%%%%%%%%%%%%%%%%%
%% collapse
load CollapseResults
DataCollapse = squeeze(mean(Results(1:SigInd,:,:),1));
h=figure(10);
axes('Parent',h,'Position',[0.1369 0.1475 0.7681 0.7775],...
    'LineWidth',1.5,...
    'FontSize',9);
image(Channel_axis,SNR_axis,256*DataCollapse')
colormap(gray(256))
colorbar('YTickLabel',{'0.2','0.4','0.6','0.8','1'})
xlabel('\# sampling channels $(\bar{m})$','Interpreter','latex','FontSize',10);
ylabel('SNR (dB)','Interpreter','latex','FontSize',10);
ax=get(h,'CurrentAxes');
set(ax,'YDir','normal')
axis(  [  4 20  -20 30 ] )


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shifts
LineWidth = 1.5;
load Results_100Shift1
DataFullRand = squeeze(mean(Results(1:SigInd,:,:),1));
load Results_20Shift5
Data20RS5 = squeeze(mean(Results(1:SigInd,:,:),1));
load Results_10Shift5
Data10RS5 = squeeze(mean(Results(1:SigInd,:,:),1));
load Results_4Shift5
Data4RS5 = squeeze(mean(Results(1:SigInd,:,:),1));
h=figure(2);
axes('Parent',h,'Position',[0.1369 0.1475 0.7681 0.7775],...
    'LineWidth',1.5,...
    'FontSize',9);
SNRInd25 = find(SNR_axis==25);
plot(Channel_axis,DataFullRand(:,SNRInd25),'b-',...
    Channel_axis,Data20RS5(:,SNRInd25),'b--',...    
    Channel_axis,Data10RS5(:,SNRInd25),'k--',...    
    Channel_axis,Data4RS5(:,SNRInd25),'r-.')   
xlabel('\# sampling channels ($\bar{m}$)','Interpreter','latex','FontSize',10);
ylabel('Empricial Recovery Rate','Interpreter','latex','FontSize',10);
axis([min(Channel_axis) max(Channel_axis) -0.05 1.1])
grid
set(h,'name','SNR = 25 dB');
ah = gca;
ach = get(ah,'Children');
for axes_i=ach
    set(axes_i,'LineWidth',LineWidth);
end

%%% -------------------------
h=figure(3);
SNRInd10 = find(SNR_axis==10);
axes('Parent',h,'Position',[0.1369 0.1475 0.7681 0.7775],...
    'LineWidth',1.5,...
    'FontSize',9);
plot(Channel_axis,DataFullRand(:,SNRInd10),'b-',...
    Channel_axis,Data20RS5(:,SNRInd10),'b--',...    
    Channel_axis,Data10RS5(:,SNRInd10),'k--',...    
    Channel_axis,Data4RS5(:,SNRInd10),'r-.')   
xlabel('\# sampling channels ($\bar{m}$)','Interpreter','latex','FontSize',10);
ylabel('Empricial Recovery Rate','Interpreter','latex','FontSize',10);
axis([min(Channel_axis) max(Channel_axis) -0.05 1.1])
grid
set(h,'name','SNR = 10 dB');
ah = gca;
ach = get(ah,'Children');
for axes_i=ach
    set(axes_i,'LineWidth',LineWidth);
end



