function [W,s]=fdla_asymmetric_sdp(P)
% P is adjacency matrix
L=size(P,1);
eye_L=eye(L);
P=P+eye_L;
zero_ind=~P;
vec_one=ones(L,1);
W_opt=ones(L,L)/L;
cvx_begin sdp quiet
    variable W(L,L)
    variable s
    minimize(s)
    subject to
        [s*eye_L,W-W_opt;
            W'-W_opt,s*eye_L]>=0;
        W*vec_one==vec_one;
        vec_one'*W ==vec_one';
        W(zero_ind)==0;
cvx_end
W(zero_ind)=0;
W=W-diag(sum(W+W')/2-ones(1,L));
end