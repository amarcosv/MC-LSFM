load '.\results\simulationData.mat'
load .\results\detectorCounts_z-0.1mm.mat
close all
NA=[0.2, 0.04, 0.025, 0.01, 0.005, 1];
NA=0.05;
image = computeNA(detectorCounts,NA,xSize/detSize,ySize/detSize);
figure,
nNA=length(NA);
for i=1:nNA
%     subplot(nNA/2,nNA/2,i)
    imagesc([xAxis(1),xAxis(end)],[yAxis(1),yAxis(end)],image(:,:,i)')
    axis equal
    xlabel('x (mm)');
    ylabel('y (mm)');
    title(sprintf('NA %.2g',NA(i)))
end

detectorCounts.ang=asin(detectorCounts.v(:,3));
detectorCounts.detid=detectorCounts.detid-1;

acceptedCounts=detectorCounts.ang>acos(NA(i)/1.33);
detectorCountsNA.detid=detectorCounts.detid(acceptedCounts);
detectorCountsNA.ppath=detectorCounts.ppath(acceptedCounts);
detectorCountsNA.v=detectorCounts.v(acceptedCounts,:);
detectorCountsNA.ang=detectorCounts.ang(acceptedCounts);
 detectorCountsNA.weights=detectorCounts.weights(acceptedCounts);

 
% detectorCounts=detectorCountsNA;

r=detectorCounts.v;

 w=detectorCounts.weights;
detPhotons=length(r);


%% Medium
% Volume size
% xSize=200; %in px
% ySize=200; %in px
% zSize=500; %in px
% voxelSize=0.1; %in mm/px

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
detectors=detMatrix(:,1:2).*dx-max(x);

detImage=zeros(xSize,ySize);

lambda=0.5*10^-6; %wavelength [m]
 lambda=1;
k=2*pi/lambda;
% k=2;


%E(r-r0)=exp[jk(a(x-x0)+b(y-y0)+c(z-z0))]

zf=0.025; %[m]
z=zf;
L1=L*dx;
for i=1:500:detPhotons
    
    r0=[detectors(detectorCounts.detid(i)-1,:),0];%%%%%%%%%%%%%%%Have to put here the exponential expression
    rX=X-r0(1);
    rY=Y-r0(2);
    rZ=X.*0;
    
     e=w(i).*exp(1j*k.*(r(i,1).*rX+r(i,2).*rY+r(i,3).*rZ)- 1j*detectorCounts.ppath(i));

%      e=exp(1j*k.*(r(i,1).*rX+r(i,2).*rY+r(i,3).*rZ));
    detImage=detImage+e;
end
I1=real(detImage); %src irradiance
%
figure(2)
imagesc(x,y,I1);
axis square; axis xy;
xlabel('x (m)'); ylabel('y (m)');
title('z= 0 m');

figure(3)
% imagesc(x,y,I2);
imagesc(x,y,real(fftshift(ifft2(ifftshift(detImage)))));
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title('I2= 2000 m');

     [focusedImage]=mcFocus(detImage,L1,lambda,zf);

% [u1]=mcFocus(detImage,L1,lambda,zf);
% figure(3)
% % imagesc(x,y,I2);
% imagesc(x,y,real(u1));
% axis square; axis xy;
% colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
% title('I2= 2000 m');

% u2=propTF(u1,L1,lambda,zf); %propagation
u2=rayProp(focusedImage,L1,lambda,zf); %propagation

I2=abs(u2.^2); %src irradiance

figure(4)
% imagesc(x,y,I2);
imagesc(x,y,real(u2));
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title(['z= ',num2str(z),' m'])

x2=x; %obs coords
y2=y;
I2=abs(u2.^2); %obs irrad

% figure(2) %display obs irrad
% imagesc(x2,y2,I2);
% axis square; axis xy;
% colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
% title(['z= ',num2str(z),' m']);
%
% figure(4) %irradiance profile
% plot(x2,I2(M/2+1,:));
% xlabel('x (m)'); ylabel('Irradiance');
% title(['z= ',num2str(z),' m']);
% %
% figure(5) %plot obs field mag
% plot(x2,abs(u2(M/2+1,:)));
% xlabel('x (m)'); ylabel('Magnitude');
% title(['z= ',num2str(z),' m']);
% %
% figure(6) %plot obs field phase
% plot(x2,unwrap(angle(u2(M/2+1,:))));
% xlabel('x (m)'); ylabel('Phase (rad)');
% title(['z= ',num2str(z),' m']);

% figure
% imagesc(x,y,log(abs(I1)));
% axis square; axis xy;
% colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
% title('z= 0 m');


% figure(3)
% imagesc(x,y,real(u2));
% axis square; axis xy;
%  xlabel('x (m)'); ylabel('y (m)');
% title('I2= 2000 m');
% 
