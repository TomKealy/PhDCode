% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

function K=getReflectionKernel(light,sensor,dim,para)
%% Compute the reflection kernel for one fixture-sensor pair
%    light: 3D spatial coordinates of the light (the fixture)
%    sensor: 3D spatial coordinates of the sensor
%    dim: 3D dimension of the room
%    para: parameter of the model
%      0: non-Lambertian
%      1: Lambertian
%    K: the resulting reflection kernel

K=zeros(dim(1),dim(2));
lx=light(1);
ly=light(2);
lz=light(3);
sx=sensor(1);
sy=sensor(2);
sz=sensor(3);

for x=1:dim(1)
    for y=1:dim(2)
        d1=sqrt((lx-x)^2+(ly-y)^2);
        d2=sqrt((sx-x)^2+(sy-y)^2);
        D1=sqrt(d1^2+lz^2);
        D2=sqrt(d2^2+sz^2);
        cos1=lz/D1;
        cos2=sz/D2;
        theta1=acos(cos1);
        v=lightDistribution(theta1)*cos1*cos2/D1^2/D2^2;
        if para==1
            v=v*cos2; % Lambertian reflectance
        end
        K(x,y)=v;
    end
end
