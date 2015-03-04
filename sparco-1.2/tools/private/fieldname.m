function str = fieldname(name)
%
% Use by: parseParams, getOption
%

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: fieldname.m 1040 2008-06-26 20:29:02Z ewout78 $

str = strrep(name,'-','_');
str = strrep(str,' ','_');
