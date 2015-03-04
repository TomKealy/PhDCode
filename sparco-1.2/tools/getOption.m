function v = getOption(opts, field, default)
%GETOPTION  Get structure field with default
%
%   GETOPTION(OPTS, FIELD, DEFAULT) returns the value of
%   OPTS.FIELD, if OPTS does not have the given FIELD, the DEFAULT
%   value is returned.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: getOption.m 1040 2008-06-26 20:29:02Z ewout78 $

  if isfield(opts,fieldname(lower(field)))
    v = getfield(opts,fieldname(lower(field)));
  else
    v = default;
  end;
