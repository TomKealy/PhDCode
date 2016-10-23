function [W,s]=fdla_asymmetric(P)
% P is adjacency matrix
L=size(P,1);
P=P+eye(L);
zero_ind=~P;
vec_one=ones(L,1);
cvx_begin quiet
    variable W(L,L)
    minimize(norm(W-ones(L,L)/L,2))
    subject to
        W*vec_one==vec_one;
        vec_one'*W ==vec_one';
        W(zero_ind)==0;
cvx_end
s=norm(W-ones(L,L)/L,2);
W(zero_ind)=0;
W=W-diag(sum(W+W')/2-ones(1,L));
end
