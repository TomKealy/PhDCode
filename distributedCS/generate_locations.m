function XY = generate_locations(number_of_locations, l);
%Generates number_of_locations points uniformly at random from a l*l
%area, points have integer co-ordinates.
XY = zeros(50,2)
XY(:,1) = randi(l,[1,number_of_locations]);
XY(:,2) = randi(l,[1,number_of_locations]);
end

