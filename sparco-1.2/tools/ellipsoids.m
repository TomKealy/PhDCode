function [P,x,y,z] = ellipsoids(varargin)
%ELLIPSOIDS  Three-dimensional Shepp-Logan phantom
%
%   [P,x,y,z] = ELLIPOIDS(G,N,OFFSET,ANGLE) creates a slice of the
%   3D Shepp-Logan phantom of size N by N. The plane is generated
%   by taking the [-1,1] x [-1,1] grid and rotating it along the z,
%   y and x axis respectively, as given by ANGLE = (X,Y,Z), and
%   recentering it to OFFSET = (X0,Y0,Z0). The angle is given in
%   radians. The ANGLE and OFFSET parameters are optional.
%
%   Parameter G can either  be a  list of ellipoids in the form 
%
%      Column 1: Intensity of the ellipse
%      Column 2: radius of the ellipsoid along the x-axis
%      Column 3: radius of the ellipsoid along the y-axis
%      Column 4: radius of the ellipsoid along the z-axis
%      Column 5: x-coordinate of the ellipsoid center
%      Column 6: y-coordinate of the ellipsoid center
%      Column 7: z-coordinate of the ellipsoid center
%      Column 8: rotation (in degrees) of the ellipse along the
%                z-axis
%
%   or one of the following four strings:
%
%      'DSLP'                    Phantom as defined in the paper by
%                                Yu, Zhao and Wang.
%
%      'Carvalho'                Phantom as defined in the PhD
%                                thesis by Carvalho
%
%      '3D Shepp-Logan'          Combination of the data given by
%                                the paper by Yu, Zhao and Wang,
%                                Carvalho's thesis and the
%                                coefficients for the 2D
%                                Shepp-Logan phantom. Taking a
%                                horizontal slice through the
%                                origin gives the standard
%                                phantom. 
%                                
%      'Modified 3D Shepp-Logan' Similar to '3D Shepp-Logan' but
%                                with increased contrast in the
%                                line with Toft. Taking a
%                                horizontal slice through the
%                                origin gives the standard modified
%                                phantom.
%
%   When omitted, G is set to 'Modified 3D Shepp-Logan'.
%
%   References
%      P. A. Toft, "The Radon Transform, Theory and Implementation",
%               PhD Thesis, 1996, (p. 199)
%      H.Y. Yu, S.Y. Zhao and G. Wang, "A Differentiable Shepp-
%               Logan Phantom and its Applications in Exact
%               Cone-beam CT", Phys. Med. Biol. 50, 2005, p.5589.
%      A. K. Jain, "Fundamentals of Digital Image Processing",
%               Prentice-Hall, NJ, 1989 (p. 439)
%      B.M. Carvalho, "Cone-beam Helict CT Virtual Endoscopy:
%               Reconstruction, Segmentation and Automatic
%               Navigation", Phd Thesis, 2003, p.50.
%
%   See also: ELLIPSES.
     
%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: ellipsoids.m 1040 2008-06-26 20:29:02Z ewout78 $

i = 1;
if (i <= nargin & ...
    (ischar(varargin{i}) || ...
     (isnumeric(varargin{i}) & size(varargin{i},2) == 8)))
  g = varargin{i}; i = i + 1;
else
  g = 'Modified 3D Shepp-Logan'; 
end

if (i <= nargin & isscalar(varargin{i}))
  n = varargin{i}; i = i + 1;
else
  n = 256;
end

if (i <= nargin & ...
    (isnumeric(varargin{i}) & length(varargin{i}) == 3))
  offset = varargin{i}; i = i + 1;
else
  offset = [0,0,-0.25];
end

if (i <= nargin & ...
    (isnumeric(varargin{i}) & length(varargin{i}) == 3))
  angle = pi * varargin{i} / 180.0;
  i = i + 1;
else
  angle = [0,0,0];
end


if ~isnumeric(g)  % z = -0.25 gives result
  switch lower(g)
   case {'DSLP'}
     % H.Y. Yu, S.Y. Zhao and G. Wang, "A Differentiable Shepp-
     % Logan Phantom and its Applications in Exact Cone-beam CT",
     % Phys. Med. Biol. 50, 2005, p.5589.
     %
     % Note: The second ellipse seems to have an incorrect y0 value
     %
     %     A     a       b      c       x0     y0      z0      phi
     g = [ 2     0.6900  0.900  0.900   0      0       0       0
          -0.98  0.6792  0.882  0.882   0      0       0       0
          -0.02  0.4100  0.160  0.210  -0.22   0      -0.25    108
          -0.02  0.3100  0.110  0.220   0.22   0      -0.25    72
           0.02  0.2100  0.250  0.500   0      0.35   -0.25    0
           0.02  0.0460  0.046  0.046   0      0.1    -0.25    0
           0.01  0.0460  0.023  0.020  -0.08  -0.65   -0.25    0
           0.01  0.0460  0.023  0.020   0.06  -0.65   -0.25    90
           0.02  0.0560  0.040  0.100   0.06  -0.105   0.625   90
          -0.02  0.0560  0.056  0.100   0      0.10    0.625   0];

     % Unused smoothness data
     g(:,9:11) = [20 40 0.98
                  20 40 0.98
                  3  6  0.90
                  3  6  0.90
                  3  6  0.90
                  2  4  0.80
                  2  4  0.80
                  2  4  0.80
                  2  4  0.80
                  2  4  0.80];
     
   case {'Carvalho'}
     % B.M. Carvalho, "Cone-beam Helict CT Virtual Endoscopy:
     % Reconstruction, Segmentation and Automatic Navigation",
     % (unpublished dissertation), 2003, p.50.
     %
     % Note: Some obvious typos were corrected.
     %
     %     A     a       b      c       x0     y0      z0      phi
     g = [ 2.0   0.6900  0.920  0.900   0      0       0       0
          -0.98  0.6624  0.874  0.880   0     -0.0184  0       0
          -0.02  0.4100  0.160  0.210  -0.22   0      -0.25    72
          -0.02  0.3100  0.110  0.220   0.22   0      -0.25   -72
           0.01  0.2100  0.250  0.350   0      0.35   -0.25    0
           0.01  0.0460  0.046  0.046   0      0.1    -0.25    0
           0.01  0.0460  0.023  0.020  -0.08  -0.605  -0.25    0
           0.01  0.0460  0.023  0.020   0.06  -0.605  -0.25    90
           0.02  0.0560  0.040  0.100   0.06  -0.105   0.625   90
          -0.02  0.0560  0.056  0.100   0      0.1     0.625   0
           0.01  0.0460  0.046  0.046   0     -0.1    -0.25    0
           0.01  0.0230  0.023  0.023   0     -0.605  -0.25    0];
           
   
   case {'3d shepp-logan'}
     % Based on
     % [1] H.Y. Yu, S.Y. Zhao and G. Wang, "A Differentiable Shepp-
     %     Logan Phantom and its Applications in Exact Cone-beam CT",
     %     Phys. Med. Biol. 50, 2005, p.5589.
     % [2] B.M. Carvalho, "Cone-beam Helict CT Virtual Endoscopy:
     %      Reconstruction, Segmentation and Automatic Navigation",
     %      (unpublished dissertation), 2003, p.50.
     % [3] A.K. Jain, "Fundamentals of Digital Image Processing",
     %     p. 439.
     %
     %     A     a       b      c       x0     y0      z0      phi
     g = [ 1     0.6900  0.920  0.900   0      0      -0.25    0
          -0.98  0.6624  0.874  0.882   0     -0.0184 -0.25    0
          -0.02  0.4100  0.160  0.210  -0.22   0      -0.25    108
          -0.02  0.3100  0.110  0.220   0.22   0      -0.25    72
           0.01  0.2100  0.250  0.500   0      0.35   -0.25    0
           0.01  0.0460  0.046  0.046   0      0.1    -0.25    0
           0.01  0.0460  0.046  0.046   0     -0.1    -0.25    0
           0.01  0.0460  0.023  0.020  -0.08  -0.605  -0.25    0
           0.01  0.0460  0.023  0.020   0.06  -0.605  -0.25    90
           0.01  0.0230  0.023  0.023   0     -0.606  -0.25    0
           0.02  0.0560  0.040  0.100   0.06  -0.105   0.625   90
          -0.02  0.0560  0.056  0.100   0      0.10    0.625   0];
     
   case {'modified 3d shepp-logan'}
     % Based on
     % [1] H.Y. Yu, S.Y. Zhao and G. Wang, "A Differentiable Shepp-
     %     Logan Phantom and its Applications in Exact Cone-beam CT",
     %     Phys. Med. Biol. 50, 2005, p.5589.
     % [2] B.M. Carvalho, "Cone-beam Helict CT Virtual Endoscopy:
     %     Reconstruction, Segmentation and Automatic Navigation",
     %     (unpublished dissertation), 2003, p.50.
     % [3] P.A. Toft, "The Random Transform, Theory and
     %     Implementation", (unpublished dissertation), 1996, p. 439.
     %
     %     A     a       b      c       x0     y0      z0      phi
     g = [ 1     0.6900  0.920  0.900   0      0      -0.25    0
          -0.8   0.6624  0.874  0.882   0     -0.0184 -0.25    0
          -0.2   0.4100  0.160  0.210  -0.22   0      -0.25    108
          -0.2   0.3100  0.110  0.220   0.22   0      -0.25    72
           0.1   0.2100  0.250  0.500   0      0.35   -0.25    0
           0.1   0.0460  0.046  0.046   0      0.1    -0.25    0
           0.1   0.0460  0.046  0.046   0     -0.1    -0.25    0
           0.1   0.0460  0.023  0.020  -0.08  -0.605  -0.25    0
           0.1   0.0460  0.023  0.020   0.06  -0.605  -0.25    90
           0.1   0.0230  0.023  0.023   0     -0.606  -0.25    0
           0.2   0.0560  0.040  0.100   0.06  -0.105   0.625   90
          -0.2   0.0560  0.056  0.100   0      0.10    0.625   0];
   
   otherwise
     g = [];
  end
end

% Compute the coordinate grid
x = repmat(linspace(-1,1,n), n,1);
y = repmat(linspace(1,-1,n)',1,n);
z = zeros(n,n);

% Rotate coordinates along z-axis
if angle(3) ~= 0
   xc = cos(angle(3)) * x + sin(angle(3)) * y;
   yc = cos(angle(3)) * y - sin(angle(3)) * x;
   x  = xc;
   y  = yc;
end

% Rotate coordinates along y-axis
if angle(2) ~= 0
  xc = cos(angle(2)) * x + sin(angle(2)) * z;
  zc = cos(angle(2)) * z - sin(angle(2)) * x;
  x  = xc;
  z  = zc;
end

% Rotate coordinates along x-axis
if angle(1) ~= 0
  zc = cos(angle(1)) * z + sin(angle(1)) * y;
  yc = cos(angle(1)) * y - sin(angle(1)) * z;
  y  = yc;
  z  = zc;
end

x = x + offset(1);
y = y + offset(2);
z = z + offset(3);

P = zeros(n,n);
for i=1:size(g,1)
   A = g(i,1);   x0 = g(i,5);
   a = g(i,2);   y0 = g(i,6);
   b = g(i,3);   z0 = g(i,7);
   c = g(i,4);   ph = g(i,8);
  
   % ------------------------------------------------------------ %
   % Equation for 3D ellipsoid centered at (x0,y0,z0) without     %
   % rotation is given by                                         %
   %                                                              %
   %          (x - x0)^2   (y - y0)^2   (z - z0)^2                %
   %          ---------- + ---------- + ---------- <= 1           %
   %             a^2          b^2          c^2                    %
   %                                                              %
   % Rotation can be added by using a rotated coordinate system.  %
   % ------------------------------------------------------------ %
   cp = cos(pi * ph / 180.0);
   sp = sin(pi * ph / 180.0);

   xc = x - x0;
   yc = y - y0;
   zc = z - z0;
   P  = P + A * ((((xc.*cp + yc.*sp).^2) / a^2 + ...
                  ((yc.*cp - xc.*sp).^2) / b^2 + ...
                  ((zc             ).^2) / c^2) <= 1);
  end
end

