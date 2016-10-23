function [W,BW,s]=fdla_symmetric_nneg(P)
% P is adjacency matrix
L=size(P,1);
P=P+eye(L);
zero_ind=~P;
vec_one=ones(L,1);
cvx_begin quiet
    variable W(L,L) symmetric
    minimize(norm(W-ones(L,L)/L,2))
    subject to
        W*vec_one==vec_one;
        W(zero_ind)==0;
        W>=0
cvx_end
s=norm(W-ones(L,L)/L,2);
W=(W+W')/2;
W(zero_ind)=0;
W=W-diag(sum(W)-ones(1,L));
BW=0.5*(eye(L)+W);
end