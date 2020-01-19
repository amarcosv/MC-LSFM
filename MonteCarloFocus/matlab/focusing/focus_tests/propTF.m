function [u2] = propTF(u1,L,lambda,z)
%PROPTF - transfer function approach
% assumes same x and y side lengths and uniform sampling
%
% Syntax:  [u2] = propTF(u1,L1,lambda,z)
%
% Inputs:
%    U1 - source plane field
%    L1 - source and observation plane side length
%    LAMBDA - wavelength
%    Z - propagation distance
%
% Outputs:
%    U2 - observation plane field
%
% Example:
%
% See also:
%
% $Author: Asier Marcos Vidal $    $Date: 11-Dec-2018$    $Revision: 0.1 $
% Copyright: 
%           BiiG - Biomedical Imaging and Instrumentation Group
%           UC3M - Universidad Carlos III de Madrid
%----------------------------- BEGIN CODE ---------------------------------

[M,N]=size(u1); %get input field array size
dx=L/M; %sample interval
k=2*pi/lambda; %wavenumber

fx=-1/(2*dx):1/L:1/(2*dx)-1/L; %freq coords
[FX,FY]=meshgrid(fx,fx);

H=exp(-j*pi*lambda*z*(FX.^2+FY.^2)); %trans func
H=fftshift(H); %shift trans func
U1=fft2(fftshift(u1)); %shift, fft src field
U2=H.*U1; %multiply
u2=ifftshift(ifft2(U2)); %inv fft, center obs field
end