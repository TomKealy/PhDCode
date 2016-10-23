%% function for model order detection 
function q= mdl_function(R,M)
% R is sample correlation matrix
% M is the number of samples that R has computed 
%% 
ev=eig(R);%eigenvalues
ev=sort(ev,'descend');
p=length(ev);%number of eigenvalues
%% detection with mdl
for m=1:p-1
    gm=geomean(ev(m+1:end));%geometrical mean 
    am=mean(ev(m+1:end));%arithmetic mean
    mdl_ic(m)=-(p-m)*M*log10(gm/am)+0.5*m*(2*p-m)*log10(M);
end    
[ic,k_min]=min(abs(mdl_ic));
q=k_min;
end
