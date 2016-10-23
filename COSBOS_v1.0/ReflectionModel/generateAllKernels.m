% Copyright (C) 2014 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

function K=generateAllKernels(lights,sensors,dim,para)
%% Generate all reflection kernels
%    lights: 3D spatial coordinates of all lights (fixtures)
%    sensors: 3D spatial coordinates of all sensors
%    dim: 3D dimension of the room
%    para: parameter of the model
%      0: non-Lambertian
%      1: Lambertian
%  Note: this function is very slow, thus we suggest you store the results
%  in a mat file

K=cell(size(sensors,1),size(lights,1));

lights(:,3)=lights(:,3);
sensors(:,3)=sensors(:,3);

for s=1:size(sensors,1)
    for l=1:size(lights,1)
        K{s,l}=getReflectionKernel(lights(l,:),sensors(s,:),dim,para);
        K{s,l}=K{s,l}(:,end:-1:1); % to be consistent with occupancy scenarios in the papers
    end
end
        