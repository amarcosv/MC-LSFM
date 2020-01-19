function [zp] = rayProp(p,v,z)
%RAYPROP - Use the parametric equations to propagate a ray in space 
%
% Syntax:  [zp] = rayProp(p,v,z)
%
% Inputs:
%    P - inital point x,y
%    V - director vector
%    Z - z end position
%
% Outputs:
%    ZP - x, y coordinates of the ray at z
%
% Example:
%
% See also:
%
% $Author: Asier Marcos Vidal $    $Date: 11-Jan-2019$    $Revision: 0.1 $
% Copyright: 
%           BiiG - Biomedical Imaging and Instrumentation Group
%           UC3M - Universidad Carlos III de Madrid
%----------------------------- BEGIN CODE ---------------------------------

% Define parametric equations
% x = x0 + t*vx;
% y = y0 + t*vy;
% z = z0 + t*vz;

t = z/v(3);
x = p(1)+ t*v(1);
y = p(2)+ t*v(2);

zp=[x,y];