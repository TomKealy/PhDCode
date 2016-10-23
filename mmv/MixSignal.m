function MixedSig = MixSignal(Sig,t_axis,Pattern,Tp)
% random mixing
M = length(Pattern);
Ind =   mod(  floor(t_axis/(Tp/M)) , M) + 1;
MixedSig = Sig .* Pattern(Ind);