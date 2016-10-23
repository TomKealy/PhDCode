%dc01e03.m
%plots the CTFT spectra of rectangular/triangular waves
clear, clf
global D
D=1; 
CTFT('rD',D,221);
CTFT('tri',D,223);
