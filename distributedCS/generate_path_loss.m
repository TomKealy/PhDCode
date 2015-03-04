function [ losses ] = generate_path_loss(distances, wavelengths)
%Uses the Friis transmission equation to generate path losses
% 
losses = zeros(50,1);
for i=1:size(losses);
    losses(i) = 20*log10((4*pi*distances(i))/wavelengths(i));
end
end

