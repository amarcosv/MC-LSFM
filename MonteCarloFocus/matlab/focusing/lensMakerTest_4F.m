clear
clc
% close all
% load results/simulation1/detectorCounts_z0mm.mat
load 'C:\Users\Asier Marcos\Downloads\testData\testData/data2fluoro.mat'
% load results/detectorCounts_z0mm.mat
% load results/simulationData.mat

%Lensmaker
MCvoxelSize=0.1;
%spimVolume=zeros(1000,1000,nzplanes);
  nzplanes=1;
zPlanesmm=0;

%photons.p(:,3)=photons.p(:,3)+100;
photons.p(:,1)=photons.p(:,1)-50;
photons.p(:,2)=photons.p(:,2)-50;
%photons.v(:,3)=-photons.v(:,3);
photons.p=photons.p.*MCvoxelSize;
%Create z plane names
for z=1:nzplanes
%zname=sprintf('detectorCounts_z%.1fmm.mat',zPlanesmm(z))  ;  
%datafile=['results/2fluoTest/' zname]

%load(datafile)

f1=100;
f2=100;
M=f2/f1;

d0=5+zPlanesmm(z);
dLens=f1-d0;

zf=f1;
f=f1;


%estimate detector coords, then referenciate to the center of the detector

x=(linspace(0,100,100)-50).*MCvoxelSize;
y=(linspace(0,100,100)-50).*MCvoxelSize;
[X,Y]=meshgrid(x,y);
r=sqrt(max(x)^2+max(y)^2);

% fprintf('Lens dimensions: %f x %f [mm]\n', maxX, maxY);
% fprintf('\t\t\t\t\t %i x %i [px]\n', detPerLine, detPerLine);
fprintf('Lens radius: \t %f mmm\n\n', r);
fprintf('Lens NA:\t\t\t%f  \n', r/f);
fprintf('Focusing at:\t\t\t%f  \n', zf);

%detPhotons=length(detPos(detectorCounts.detid));
% detectorCounts.v=cos(detectorCounts.v);
%   detPhotons=10000000;
thetaDev=sqrt(X.^2+Y.^2)./f;
rDev=f-sqrt(X.^2+Y.^2);
rDev=rDev-max(max(rDev));
rDev=rDev.*0;


%Create a plane for the detector
%Define our camera detector
pxXd=200;
pxYd=200;
detPlane=zeros(pxXd,pxYd);

detPlane2=zeros(pxXd,pxYd);
pxSized=0.05; %in mm/px
xd=0:pxXd-1; %x coords
yd=0:pxYd-1; %y coords
xd=xd.*pxSized-pxXd*pxSized/2;
yd=yd.*pxSized-pxYd*pxSized/2;
[Xd,Yd]=meshgrid(xd,yd); %[in mm units]


focusedV=zeros(Nphotons,3);
detPlanePhotons=zeros(Nphotons,2);
% detPlanePhotons2=zeros(Nphotons,2);
% detPlane=zeros(pxX,pxY);

% x=0:pxX-1; %x coords
% y=0:pxY-1; %y coords
% x=x.*pxSize;
% y=y.*pxSize;
% [X,Y]=meshgrid(x,y);

fprintf('Detector dimensions: %f x %f [mm]\n', pxXd*pxSized, pxYd*pxSized);
fprintf('\t\t\t\t\t %i x %i [px]\n', pxXd, pxYd);

% figure(1)
% xlabel('mm');
% ylabel('mm');
% figure(2)
% xlabel('mm');
% ylabel('mm');
Nphotons=100000;
for i=1:1:Nphotons
    lensPos=rayProp(photons.p(i,:),photons.v(i,:),dLens);
    
    focusedV(i,:)=rotateVect(lensPos,photons.v(i,:),-sqrt(lensPos(1).^2+lensPos(2).^2)./f1);
    lensPos=rayProp(lensPos,focusedV(i,:),f1+f2);
    focusedV(i,:)=rotateVect(lensPos, focusedV(i,:),-sqrt(lensPos(1).^2+lensPos(2).^2)./f2);
    
    detPlanePhotons(i,:)=rayProp(lensPos,focusedV(i,:),f2);
    
    detPlanePhotons(i,:)=round(detPlanePhotons(i,:)./pxSized)+pxXd/2+1;
    if detPlanePhotons(i,:)<pxXd & detPlanePhotons(i,:)>0
        detPlane(detPlanePhotons(i,1),detPlanePhotons(i,2))=detPlane(detPlanePhotons(i,1),detPlanePhotons(i,2))+photons.w(i);
        
    end
    
end
detPlanePhotons=zeros(Nphotons,2);
% for i=1:1:Nphotons
%     lensPos=rayProp(photons.p(i,:),photons.v(i,:),dLens);
%     lensDeltaf=abs(sqrt(f1^2-(lensPos(1).^2+lensPos(2)))-f1);
%      lensPos=rayProp(lensPos,photons.v(i,:),-lensDeltaf);
% %     focusedV(i,:)=rotateVect(lensPos,photons.v(i,:),-sqrt(lensPos(1).^2+lensPos(2).^2)./f1);
%      focusedV(i,:)=rotateVect(lensPos,photons.v(i,:),-thetaCalc(sqrt(lensPos(1).^2+lensPos(2).^2),f1));   
%     lensPos=rayProp(lensPos,focusedV(i,:),f1+f2+lensDeltaf);
%     lensDeltaf=abs(sqrt(f2^2-(lensPos(1).^2+lensPos(2)))-f2);
%     lensPos=rayProp(lensPos,focusedV(i,:),lensDeltaf);
% %     focusedV(i,:)=rotateVect(lensPos, focusedV(i,:),-sqrt(lensPos(1).^2+lensPos(2).^2)./f2);
%     focusedV(i,:)=rotateVect(lensPos, focusedV(i,:),-thetaCalc(sqrt(lensPos(1).^2+lensPos(2).^2),f2));
%     
%     detPlanePhotons(i,:)=rayProp(lensPos,focusedV(i,:),f2-lensDeltaf);
%     
%     detPlanePhotons(i,:)=round(detPlanePhotons(i,:)./pxSized)+pxXd/2+1;
%     if detPlanePhotons(i,:)<pxXd & detPlanePhotons(i,:)>0
%         detPlane2(detPlanePhotons(i,1),detPlanePhotons(i,2))=detPlane2(detPlanePhotons(i,1),detPlanePhotons(i,2))+photons.w(i);        
%     end    
% end

%spimVolume(:,:,z)=detPlane;

% figure,
% imagesc(xd./M,yd./M,detPlane);
% colorbar
% 
% figure,
% imagesc(xd./M,yd./M,log(detPlane+1));
% colorbar


detPlane=detPlane./max(max(detPlane));
% detPlane2=detPlane2./max(max(detPlane2));
%% Now make a beautiful plot
figure,
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(gcf,'color','w');

subplot(1,2,1)
imagesc(xd./M,yd./M,detPlane);
xlabel('x [mm]','Interpreter','latex');
ylabel('y [mm]','Interpreter','latex');
colorbar
% cmocean('haline')
pbaspect([1 1 1])
title('Fluorescence image','Interpreter','latex')

% subplot(1,2,2)
% imagesc(xd./M,yd./M,detPlane2);
% xlabel('x [mm]','Interpreter','latex');
% ylabel('y [mm]','Interpreter','latex');

subplot(1,2,2)
imagesc(xd./M,yd./M,log10(detPlane+1));
xlabel('x [mm]','Interpreter','latex');
ylabel('y [mm]','Interpreter','latex');
colorbar
% cmocean('haline')
 pbaspect([1 1 1])
title('Fluorescence log image','Interpreter','latex')
set(gcf,'position',[10,10,1500,750])
mua=0.5;
mus=10;
tit=sprintf('Monte Carlo SPIM Simulator [z=%.1f mm - \\mu_a=%6.2f cm^{-1}, \\mu''_s=%6.2f cm^{-1}]\n',zPlanesmm(z),mua,mus);
%suptitle(tit);
drawnow
end
% 
% 
% fileID = fopen('SPIM_vol.img','w');
% fwrite(fileID,single(spimVolume),'single');
% fclose(fileID);

