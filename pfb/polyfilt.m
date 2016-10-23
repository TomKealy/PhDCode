%
% Polyphase filter implementation (2 channels)
% 
%   X = input signal, separated into even and odd phases.
%       first row = even phase
%	second row = odd phase	
%   Y = output signal, separated into even and odd phases.
%   H = 2x2 polyphase matrix
%       H(1,1,:) = h0,even[n]
%       H(1,2,:) = h0,odd[n]
%       H(2,1,:) = h1,even[n]
%       H(2,2,:) = h1,odd[n]
%
function Y = polyfilt(H, X)
y0 = conv(H(1,1,:),X(1,:)) + conv(H(1,2,:),X(2,:));
y1 = conv(H(2,1,:),X(1,:)) + conv(H(2,2,:),X(2,:));
Y = [y0; y1];
