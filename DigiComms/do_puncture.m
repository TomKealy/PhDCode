%do_puncture.m
clear
punction_vector=[1 1 0 1 0 1]; 
integer_seq=[1:12]; % An integer sequence to puncture
punctured=puncture(integer_seq,punction_vector) % To see the punctured result
depunctured=depuncture(punctured,punction_vector);
sim('do_puncture_sim',0.99)
[integer_seq; depunctured; depunctured_sim.'] % To see the depunctured result
