close all
clear all
% close all
%% Medium
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
% x=-L/2:dx:L/2-dx; %x coords
% y=-M/2:dy:M/2-dy; %y coords
x=-L/2:1:L/2-1; %x coords
y=-M/2:1:M/2-1; %y coords
x=x.*dx;
y=y.*dy;
[X,Y]=meshgrid(x,y);
% detectors=detMatrix(:,1:2).*dx-max(x);

detImage=zeros(xSize,ySize);

  lambda=0.5*10^-6; %wavelength [m]
%   lambda =0.01;
k=2*pi/lambda;
% k=100;

% figure,
% ang=deg2rad(0:1:360);
% ang=deg2rad(90);

% imagesc(x,y,(detImage));
% axis square; axis xy;
% colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
% title('z= 0 m');
% c=1;


Etotal=zeros(xSize,ySize);
Eprop=zeros(xSize,ySize);
%% Generate here a number of random direction photons
nPhotons=1000;
phi=2*pi.*rand([1,nPhotons]);
theta=pi/2*rand([1,nPhotons]);
kx=k*cos(phi);
ky=k*sin(phi);
% kz=sqrt(k.^2-kx.^2-ky.^2);
kz=k*cos(theta);

Sx=0;
Sy=0;

Z=0;
figure,
for ph=1:nPhotons
eX=exp(1i*kx(ph).*X);
eY=exp(1j*ky(ph).*Y);
eZ=exp(1j*kz(ph).*Z);
E=eX.*eY.*eZ.*exp(1i*kx(ph).*Sx).*exp(1i*ky(ph).*Sy);

imagesc(x,y,real(E));
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
Etotal=E+Etotal;
end
figure,
imagesc(x,y,real(Etotal));
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');



Z=0.01;
for ph=1:nPhotons
[S]=rayProp([Sx,Sy,0],[kx(ph)/k,ky(ph)/k,kz(ph)/k],Z);
eX=exp(1i*kx(ph).*X);
eY=exp(1j*ky(ph).*Y);
eZ=exp(1j*kz(ph).*Z);
E=eX.*eY.*eZ.*exp(1i*kx(ph).*S(1)).*exp(1i*ky(ph).*S(2));
Eprop=E+Eprop;
end
figure,
imagesc(x,y,abs(real(Eprop)));
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');

% 
% Z=5;
% 
% 
% 
% 
% figure,
% imagesc(x,y,real(Eprop));
% axis square; axis xy;
% colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
