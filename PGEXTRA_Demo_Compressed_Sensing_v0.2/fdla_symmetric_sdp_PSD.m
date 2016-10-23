function [W,s]=fdla_symmetric_sdp_PSD(P)
% P is adjacency matrix
L=size(P,1);
eye_L=eye(L);
P=P+eye_L;
zero_ind=~P;
vec_one=ones(L,1);
W_opt=ones(L,L)/L;
cvx_begin sdp quiet
    variable W(L,L) symmetric
    variable s
    minimize(s)
    subject to
        [s*eye_L,W-W_opt;
            W-W_opt,s*eye_L]>=0;
        W*vec_one==vec_one;
        W(zero_ind)==0;
        W>=0;
cvx_end
W=(W+W')/2;
W(zero_ind)=0;
W=W-diag(sum(W)-ones(1,L));
end