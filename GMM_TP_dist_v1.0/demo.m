
clc
clear
close all

s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s)

filename = 'datasets/TrafficT8_P8X8X8/Traffic';


load(filename);
[Row,Col,nF] = size(Xtst);


[Row,Col,nF] = size(Xtst);
para.filename = filename;
sigma = 0.001;
T = 8;  % How many frames collapsed to 1 measurement
M = 12;   % How many measurements used for reconstruct

Xtst = Xtst(:,:,1:T*M);
C = binornd(1,0.5,[Row, Col,T]);
shift = 1;
for t=2:T
    C(:,1+(t-1)*shift:Col,t) = C(:,1+(t-2)*shift:Col-shift,t-1);
end
Y = zeros(Row,Col,M);
for t = 1:M
    Y(:,:,t) = sum(Xtst(:,:,(t-1)*T+(1:T)).*C,3) + sigma*ones(Row,Col);
end
save('data_global','Y', 'C','Xtst');

Algorithm = 1;
switch Algorithm
    case 1 % learn the GMM offline
        para.method = 'GMM';
        addpath('GMM/')
        Xpre = interface_GMM(para);
        rmpath('GMM/')
    case 2 % learn the GMM offline and update the GMM online. Caution: it takes longer time.
        para.method = 'GMM_online';
        addpath('GMM/')
        para.dirc = 'data0_P8X8X8';  
        Xpre = interface_GMM_online(para);
        rmpath('GMM/')
end
time = tic;
[PSNR, SSIM] = saveResults(Xpre, Row, Col, M, T,Y,Xtst,filename,para.method,time);
time = toc;
