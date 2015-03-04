function str = opToString(op)
%OPTOSTRING  Generate string description of operator
%
%   OPTOSTRING(OP) returns a string describing operator OP.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opToString.m 1040 2008-06-26 20:29:02Z ewout78 $

opInfo = op([],0);
type = opInfo{4};
switch type{1}
   case {'Dictionary','BlockDiag','Stack'}
      str = [type{1} '('];
      oplist = type{3};
      for i=1:length(oplist)
         s2 = opToString(oplist{i});
         if i ~= 1,
           str = [str, ', ', s2];
         else
           str = [str, s2];
         end
      end
      str = [str, ')'];

   case {'FoG'}
      oplist = type{2}; str = '';
      for i=1:length(oplist)
         s2 = opToString(oplist{i});
         if i~= 1,
            str = [str, ' * ', s2];
         else
            str = [str, s2];
         end
      end

   case {'Sum'}
      oplist = type{2}; str = '(';
      for i=1:length(oplist)
         s2 = opToString(oplist{i});
         if i~= 1,
            str = [str, ' + ', s2];
         else
            str = [str, s2];
         end
      end
      str = [str,')'];
 
   case {'Transpose'}
      oplist = type{2};
      str = ['[',opToString(oplist), ']^T'];

   case {'Crop'}
      str = sprintf('Crop(%s)',opToString(type{2}));
 
   case {'Wavelet'}
      name = type{2};
      name(1) = upper(name(1));
      str = sprintf('%s-Wavelet',name);
      
   case {'Kron'}
      str = sprintf('Kron(%s,%s)',opToString(type{2}), ...
                                  opToString(type{3}));
   
   case {'WindowedOp'}
      str = sprintf('Windowed(%s)',opToString(type{2}));
    
 
 otherwise
      str = type{1};
end
