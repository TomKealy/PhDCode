function x=tri(t)
global D
tmp=abs(t); x=(tmp<=D).*(1-tmp/D);
