function [figtype, figext] = getFigureExt(figtype)
%GETFIGUREEXT  Get default filename extension for print device
%
%   [FIGTYPE,FIGEXT] = GETFIGUREEXT(FIGTYPE) returns the figure
%   extension FIGEXT corresponding to FIGTYPE. For example, passing
%   device type 'epsc2' gives a FIGEXT equal to 'eps'. When FIGTYPE
%   is not recognized FIGTYPE and FIGEXT are set to 'png' by
%   default. 

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: getFigureExt.m 1040 2008-06-26 20:29:02Z ewout78 $

switch lower(figtype)
 case {'pdf'}
   figext = 'pdf';
   
 case {'png'}
   figext = 'png';
   
 case {'eps','epsc','eps2','epsc2'}
   figext = 'eps';
   
 case {'ps','psc','ps2','psc2'}
   figext = 'ps';
   
 otherwise
   figtype = 'png';
   figext  = 'png';
end
