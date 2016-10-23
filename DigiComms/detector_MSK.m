function d=detector_MSK(dth)
if dth>0, d=1; else d=0; end
if abs(dth)>pi, d=1-d; end
