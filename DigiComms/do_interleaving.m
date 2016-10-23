%do_interleaving.m
clear
%Nrow=16, Ndbps=144, code_rate=3/4, Nbpsc=4;
Nrow=8, Ndbps=96, code_rate=1/2, Nbpsc=4
Ncbps=Ndbps/code_rate, Ncol=Ncbps/Nrow; s=max(Nbpsc/2,1); 
integer_seq=[1:Ncbps]; % A pseudo coded sequence to interleave
Integer_seq=reshape(integer_seq,Nrow,Ncol) % To see the original sequence
interleaved=interleaving(integer_seq,Nrow,Ncbps,Nbpsc,1);
Interleaved=reshape(interleaved,Nrow,Ncol) % To see the interleaved result
deinterleaved=deinterleaving(interleaved,Nrow,Ncbps,Nbpsc);
Deinterleaved=reshape(deinterleaved,Nrow,Ncol) % To see deinterleaved result
sim('do_interleaving_sim',0.99) % Run the Simulink model
[interleaved; interleaved_sim] % Do the two interleaved results conform?
sum(deinterleaved~=deinterleaved_sim) % Do the deinterleaved results conform?
sum(integer_seq~=deinterleaved) % deinterleaved conforms with the original data
