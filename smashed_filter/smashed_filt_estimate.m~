function [estimate] = smashed_filt_estimate(y, M, K, A, L, F)

%
% y -> The compressive measurement
% M -> Size of original signal
% K -> Size of compressive signal
% L -> 2nd difference matrix
% F -> Lehmer matrix of size M
%

direct_templates = zeros(M, K);
c = zeros(1, M);

for k=1:M
    direct_templates(k, :) = A*L(k, :)';
    c(k) = K*(y'*direct_templates(k,:)');
end

ghat = L\c';
ahat = F\c';

% Find k biggest values of ahat

[sorted_ahat, inds] = sort(abs(ahat), 'descend');

num_change_points = 80;

ahat_big = inds(1:num_change_points);

for ii = 1:M
    if ismember(ii, ahat_big);
        ii;
    else
        ahat(ii) = 0;
    end
end

% take the average of ghat between these vaules
piece_mean = zeros(1, num_change_points);
piece_var = zeros(1,num_change_points);
threshths = zeros(1,num_change_points);
inds = [1, ahat_big', 300];
inds = sort(inds);
piece_ttests = zeros(1, num_change_points);

for l=1:num_change_points
    piece = ghat(inds(l):inds(l+1))';
    piece_mean(l) = mean(piece);
    piece_var(l) = var(piece,1);
    threshths(l) = (1.96*piece_var(l))/sqrt(length(piece));
    piece_ttests(l) = ttest(piece, 0, 'Tail', 'right');
end

estimate = zeros(1, M);

for ii=1:num_change_points
    for jj = inds(ii):inds(ii+1)
        if piece_ttests(ii) == 1
            estimate(jj) = 1; %piece_mean(ii);
        else
            estimate(jj) = 0;
        end
    end
end

end

