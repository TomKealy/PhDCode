function P = ellipses(g,n)
%ELLIPSES  Two-dimensional Shepp-Logan phantom
%
%   P = ELLIPSES(G,N) creates a Shepp-Logan phantom of size N by
%   N. Parameter G can either be a list of ellipses in the form
%
%      Column 1: Intensity of the ellipse
%      Column 2: width of the ellipse
%      Column 3: height of the ellipse
%      Column 4: x-coordinate of the ellipse center
%      Column 5: y-coordinate of the ellipse center
%      Column 6: rotation of the ellipse
%
%   or one of the following two strings:
%
%      'Shepp-Logan'           Original phantom
%      'Modified Shepp-Logan'  Phantom with improved contrast.
%
%   References
%      A. K. Jain, "Fundamentals of Digital Image Processing",
%                  Prentice-Hall, NJ, 1989 (p. 439)
%      P. A. Toft, "The Radon Transform, Theory and Implementation",
%                  PhD Thesis, 1996, (p. 199)
%
%   See also: ELLIPSOIDS.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: ellipses.m 1040 2008-06-26 20:29:02Z ewout78 $


if nargin < 2, n = 256; end;
if nargin == 1 & isscalar(g)
   n = g; g = 'modified Shepp-Logan';
end
if nargin < 1, g = 'modified Shepp-Logan'; end

P = zeros(n,n);

if ~isnumeric(g)
  switch lower(g)
   case {'modified shepp-logan'}
     %      A    a       b     x0    y0     phi
     g = [ 1.0, 0.69,   0.92,  0.0,  0.0,    0; ...
          -0.8, 0.6624, 0.874, 0.0, -0.0184, 0; ...
          -0.2, 0.11,   0.31,  0.22, 0.0,   -pi/10; ...
          -0.2, 0.16,   0.41, -0.22, 0.0,    pi/10; ...
           0.1, 0.21,   0.25,  0.0,  0.35,   0; ...  
           0.1, 0.046,  0.046, 0.0,  0.1,    0; ...
           0.1, 0.046,  0.046, 0.0, -0.1,    0; ...
           0.1, 0.046,  0.023,-0.08,-0.605,  0; ...
           0.1, 0.023,  0.023, 0.0, -0.605,  0; ...
           0.1, 0.023,  0.046, 0.06,-0.605,  0];
   case {'shepp-logan'}
     %      A    a       b     x0    y0     phi
     g = [ 1.00, 0.69,   0.92,  0.0,  0.0,    0; ...
          -0.98, 0.6624, 0.874, 0.0, -0.0184, 0; ...
          -0.02, 0.11,   0.31,  0.22, 0.0,   -pi/10; ...
          -0.02, 0.16,   0.41, -0.22, 0.0,    pi/10; ...
           0.01, 0.21,   0.25,  0.0,  0.35,   0; ...  
           0.01, 0.046,  0.046, 0.0,  0.1,    0; ...
           0.01, 0.046,  0.046, 0.0, -0.1,    0; ...
           0.01, 0.046,  0.023,-0.08,-0.605,  0; ...
           0.01, 0.023,  0.023, 0.0, -0.605,  0; ...
           0.01, 0.023,  0.046, 0.06,-0.605,  0];
    
   otherwise
     g = [];
  end
end


x = repmat(linspace(-1,1,n), n,1);
y = repmat(linspace(1,-1,n)',1,n);

for i=1:size(g,1)
  A = g(i,1);   x0 = g(i,4);
  a = g(i,2);   y0 = g(i,5);
  b = g(i,3);   ph = g(i,6);
  
  % ------------------------------------------------------------ %
  % Equation for ellipse centered at (x0,y0) without rotation is %
  % given by                                                     %
  %               (x - x0)^2   (y - y0)^2                        %
  %               ---------- + ---------- <= 1                   %
  %                  a^2          b^2                            %
  %                                                              %
  % Rotation can be added by using a rotated coordinate system.  %
  % This is done by multiplying by the rotation matrix           %
  %                                                              %
  %           |x'|     |  cos(phi)   sin(phi) | |x|              %
  %           |y'| =   | -sin(phi)   cos(phi) | |y|.             %
  %                                                              %
  % ------------------------------------------------------------ %
  cp = cos(ph);
  sp = sin(ph);
  xc = x - x0;
  yc = y - y0;
  P  = P + A * ((((xc.*cp + yc.*sp).^2) / a^2 + ...
                 ((yc.*cp - xc.*sp).^2) / b^2 ) <= 1);
end
