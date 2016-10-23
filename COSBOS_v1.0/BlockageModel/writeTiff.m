% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

function writeTiff(V,filename)
%% Write a 3D volume to a TIFF image

maxValue=max(V(:));
minValue=min(V(:));
V=(V-minValue)/(maxValue-minValue)*255;
V=round(V);
V=uint8(V);

for i=1:size(V,3)
    if i==1
        imwrite(V(:,:,i),filename);
    else
        imwrite(V(:,:,i),filename,'WriteMode','append');
    end
end