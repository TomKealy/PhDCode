function [gM,m]=gm2gM(gm)
% Convert into PNG polynomial description used by Simulink -Eq.(10.1.2)
if max(gm)>1, m=max(gm); gM=[m m-gm]; else m=length(gm); gM=[1 gm]; end