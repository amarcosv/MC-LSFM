 close all
clear all
% close all
%% Medium
% Volume size
xSize=200; %in px
ySize=200; %in px
zSize=500; %in px
voxelSize=0.011; %in mm/px

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
 lambda =0.01;
k=2*pi/lambda;
% k=100;

figure,
ang=deg2rad(0:1:360);
% ang=deg2rad(90);

imagesc(x,y,(detImage));
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title('z= 0 m');
c=1;
for i=1:length(ang)

Z=0;
phi=2*pi*rand(1);
theta=pi*rand(1);
% alpha=ang(i);
% beta=ang(length(ang)-i+1);
% alpha=sin(ang(i))*c;
% beta=cos(ang(i))*c;
% gamma=sqrt(1-alpha^2-beta^2);
%k=kx+ky+kz=k(alpha+beta+gamma);
%E(r-r0)=exp[jk(a(x-x0)+b(y-y0)+c(z-z0))]
eX=exp(1i*k*cos(phi).*(X));
eY=exp(1j*k*sin(phi).*(Y));
eZ=exp(1j*k*cos(theta).*Z);
et=exp(1j*rand(1)*pi*0.8);

 E=eX.*eY.*eZ.*et;
% E=eX.*eY.*eZ;
detImage=detImage+E;
% I1=real(E); %src irradiance
% abs
%
Et = real(detImage.*(detImage'));
I2=abs(Et.^2); %src irradiance
imagesc(x,y,real(detImage));
 imagesc(x,y,real(E));
pause(0.01)

end
figure,
% imagesc(x,y,I2);
% imagesc(x,y,abs(fftshift(fft2(detImage))));
imagesc(x,y,abs(real(detImage)));
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title('I2= 2000 m');
z=10;
% [u2] = rayProp(detImage,L*dx,lambda,z)
% figure,
% I2=abs(u2.^2); %src irradiance
% imagesc(x,y,real(I2));
% axis square; axis xy;
% colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
% title('z= 0 m');

% figure, imshow(abs(angle(eX)),[])