close all
clear all
% close all
% Medium
% Volume size
xSize=200; %in px
ySize=200; %in px
zSize=500; %in px
voxelSize=0.1; %in mm/px

% Coordinates space
% [M,N]=size(uin); %get input field array size
% dx=L/M; %sample interval
% k=2*pi/lambda; %wavenumber

L=xSize;
M=ySize;
dx=voxelSize;
dy=dx;

x=-L/2:1:L/2-1; %x coords
y=-M/2:1:M/2-1; %y coords
x=x.*dx;
y=y.*dy;
[X,Y]=meshgrid(x,y);
% detectors=detMatrix(:,1:2).*dx-max(x);

detImage=ones(xSize,ySize);



% lambda=0.5*10^-6; %wavelength [m]
% lambda =0.01;
% k=2*pi/lambda;
% k=100;
