function [PSNR,SSIM] = saveResults(Xpre, Row, Col, M, T,y,Xtst,filename,Algorim,time)
Xrecon = zeros(Row,Col,M*T);

for k=1:(M*T)
    Xrecon(:,:,k) = Xpre(:,:,k);
    PSNR(k) = psnr(Xrecon(:,:,k), Xtst(:,:,k));
    SSIM(k) = ssim(Xrecon(:,:,k), Xtst(:,:,k));
end

savename = [filename '_T' num2str(T) '_F' num2str(M) '_' Algorim];
format short
save(savename, '-v7.3', 'Xrecon','Xpre','PSNR','time','SSIM');

writerObj = VideoWriter([savename '.mp4'],'MPEG-4');
writerObj.FrameRate = 12;
open(writerObj);
scrsz = get(0,'ScreenSize');
fig=figure('Position',[50 100 floor(scrsz(3)*0.8) floor(scrsz(4)*0.6)]);
for nF=1:(M*T)
    subplot(1,3,1);
    imshow(Xtst(:,:,nF));
    title(['Original Video, Frame:' num2str(nF)]);
    
    
    subplot(1,3,2);
    imshow(y(:,:,ceil(nF/T))/max(max(y(:,:,ceil(nF/T)))));
    title(['Measurement video, Frame: ' num2str(ceil(nF/T))]);
    
    subplot(1,3,3);
    imshow(Xrecon(:,:,nF));
    title(['Recon video, PSNR: ' num2str(PSNR(nF))]);
        
    %pause(0.02);
    frame = getframe(fig);
    writeVideo(writerObj,frame);
end
close(writerObj);

