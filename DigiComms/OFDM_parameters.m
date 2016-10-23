function [CR,DR,Nbpsc,Ncbps,Ndbps,pun_vec,Mod,pps]=OFDM_parameters(RATE)
s='RATE?(Type 0/1/2/3/4/5/6/7/8' for 6(no coding)/6/9/12/18/24/36/48/54Mbps)';
if nargin<1, RATE=input(s); end
switch RATE % Fig. 11.18 and Table 11.2 - Rate-dependent parameters
  case 0, CR='1/1'; DR=6; Nbpsc=1; Ndbps=48; Mod='PSK';
  case 1, CR='1/2'; DR=6; Nbpsc=1; Ndbps=24; Mod='PSK';
  case 2, CR='3/4'; DR=9; Nbpsc=1; Ndbps=36; Mod='PSK';
  case 3, CR='1/2'; DR=12; Nbpsc=2; Ndbps=48; Mod='PSK';
  case 4, CR='3/4'; DR=18; Nbpsc=2; Ndbps=72; Mod='PSK';
  case 5, CR='1/2'; DR=24; Nbpsc=4; Ndbps=96; Mod='QAM';
  case 6, CR='3/4'; DR=36; Nbpsc=4; Ndbps=144; Mod='QAM';
  case 7, CR='2/3'; DR=48; Nbpsc=6; Ndbps=192; Mod='QAM';
  otherwise, CR='3/4'; DR=54; Nbpsc=6; Ndbps=216; Mod='QAM';
end
switch CR
  case '1/2', pun_vec=[]; % coding with gain K/N=1/2
  case '2/3', pun_vec=[1 1 1 0]; 
     % (1/2)*(4/3)=2/3 by coding (1/2) and puncturing (4/3)
  case '3/4', pun_vec=[1 1 0 1 0 1]; 
     % (1/2)*(6/4)=3/4 coding (1/2) and puncturing (6/4)
  otherwise  pun_vec=[]; 
end
K=1; N=2; % # of input/outputs of the convolutional encoder Gc=[171 133];
fprintf('\n Data_rate = %d[Mbps] with modulation = 2^%d-%s',DR,Nbpsc,Mod)
form='\n Code_rate=Original_Code_Rate/Puncture_Rate=(%d/%d)/(%d/%d)=%s';
fprintf(form, K,N,sum(pun_vec),length(pun_vec),CR);
form='\n Ncbps(# of Code Bits Per OFDM Symbol)=Ndbps/Code_Rate=%d/(%s) = %d\n';
Ncbps=Ndbps/str2num(CR);  fprintf(form, Ndbps,CR,Ncbps);
scrambler_poly = [0 0 0 1 0 0 1];  x0=ones(1,7); 
pps = rotate_l(1-2*png(scrambler_poly,0,x0),7); % Pilot polarity sequence