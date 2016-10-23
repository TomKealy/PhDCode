data = csvread('tvws_data1.csv');

data_length = 10000;
x = data(1:data_length,2);
signal = x;
A = randn(data_length, data_length);

y = A*x;

lambda = sqrt(12*2*log(data_length));
rho = nthroot(1/max(abs(eig(A'*A))), 3);

xhat = lasso(A, y, lambda, rho, 1.0);
%xbar = total_variation(y, lambda, rho, 1.0);

% %wavelet stuff
% 
% lev   = 20;
% wname = 'db1';
% nbcol = 64;
% [c,l] = wavedec(signal,lev,wname);
% 
% %Expand discrete wavelet coefficients for plot.
% 
% len = length(signal);
% cfd = zeros(lev,len);
% for k = 1:lev
%     d = detcoef(c,l,k);
%     d = d(:)';
%     d = d(ones(1,2^k),:);
%     cfd(k,:) = wkeep1(d(:)',len);
% end
% cfd =  cfd(:);
% I = find(abs(cfd)<sqrt(eps));
% cfd(I) = zeros(size(I));
% cfd    = reshape(cfd,lev,len);
% cfd = wcodemat(cfd,nbcol,'row');
% 
% %Perform the continuous wavelet transform (CWT) and visualize results
% 
% h311 = subplot(3,1,1);
% h311.XTick = [];
% plot(signal,'r');
% title('Analyzed signal.');
% ax = gca;
% ax.XLim = [1 length(signal)];
% subplot(3,1,2);
% colormap(cool(128));
% image(cfd);
% tics = 1:lev;
% labs = int2str(tics');
% ax = gca;
% ax.YTickLabelMode = 'manual';
% ax.YDir = 'normal';
% ax.Box = 'On';
% ax.YTick = tics;
% ax.YTickLabel = labs;
% title('Discrete Transform, absolute coefficients.');
% ylabel('Level');
% h312 = subplot(3,1,2);
% h312.XTick = [];
% subplot(3,1,3);
% scales = (1:32);
% cwt(signal,scales,wname,'plot');
% colormap(cool(128));
% ax = gca;
% tt = ax.YTickLabel;
% [r,c] = size(tt);
% yl = char(32*ones(r,c));
% for k = 1:3:r
%     yl(k,:) = tt(k,:);
% end
% ax.YTickLabel = yl;