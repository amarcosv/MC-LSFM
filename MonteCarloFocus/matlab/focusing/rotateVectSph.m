function [vr] = rotateVectSph(p,v,theta)
%ROTATEVECT - Fuction to rotate a vector located at P towards the z axis.
%The vector will perform a convergent or divergent rotation

%
% Syntax:  [] = rotateVect(P,v,theta)
%
% Inputs:
%    P - Position of the vector to be rotated
%    V - vector
%    THETA - rotation angle
%
% Outputs:
%    VR - rotated vector
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

% Normalice the new base
pMod=sqrt(p(1)^2+p(2)^2);
p=p./pMod;

%Create the Base Change Matrix (from P to canonical) and its inverse
MBpBc=[p(1) -p(2) 0; ...
       p(2)  p(1) 0; ...
        0     0   1 ];
%MBcBp=inv(MBpBc); Instead of this, we use the operator \

%Create rotation matrix around y (now beta) axis
Rb=[cos(theta) 0 sin(theta); ...
       0       1    0;       ...
   -sin(theta) 0 cos(theta) ];



%Change base
%vBp=MBcBp*v';
vBp=MBpBc\v';

mvBp=norm(vBp);
thetaV=acos(vBp(3)/mvBp);
thetaV=thetaV+theta;
%Rotate
% vrBp=Rb*vBp;

vrBp=[mvBp*sin(thetaV) mvBp*sin(thetaV)*sin(atan(vBp(3)/vBp(2))) mvBp*cos(thetaV)]';
%return to canonical base

vr=MBpBc*vrBp;
vrMod=sqrt(vr(1)^2+vr(2)^2+vr(3)^2);
vr=vr./vrMod;
    