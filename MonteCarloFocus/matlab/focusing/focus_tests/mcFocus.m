function [uout] = mcFocus(uin,L,lambda,zf)
%MCFOCUS - One line description of what the function or script performs (H1 line)
%Additional file header info
%
% Syntax:  [] = mcFocus()
% uniform sampling assumed
%
% Inputs:
%   uin - input field
%   L - side length
%   lambda - wavelength
%   zf - focal distance (+ converge, - diverge)
% Outputs:
%   uout - output field
% Example:
%
% See also:
%
% $Author: Asier Marcos Vidal $    $Date: 10-Dec-2018$    $Revision: 0.1 $
% Copyright: 
%           BiiG - Biomedical Imaging and Instrumentation Group
%           UC3M - Universidad Carlos III de Madrid
%----------------------------- BEGIN CODE ---------------------------------

% converging or diverging phase-front



[M,N]=size(uin); %get input field array size
dx=L/M; %sample interval
k=2*pi/lambda; %wavenumber

 x=-L/2:dx:L/2-dx; %coords
% x=-L/2:1:L/2-1; %x coords
% x=x.*dx;
[X,Y]=meshgrid(x,x);

uout=uin.*exp(-1j*k/(2*zf)*(X.^2+Y.^2)); %apply focus
 end