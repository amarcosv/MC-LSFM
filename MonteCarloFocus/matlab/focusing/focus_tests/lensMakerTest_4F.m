clear
clc
% close all
% load results/simulation1/detectorCounts_z0mm.mat
load results/simulation1/simulationData.mat
% load results/detectorCounts_z0mm.mat
% load results/simulationData.mat

%Lensmaker

spimVolume=zeros(1000,1000,nzplanes);
  nzplanes=1;
zPlanesmm=0;
%Create z plane names
for z=1:nzplanes
zname=sprintf('detectorCounts_z%.1fmm.mat',zPlanesmm(z))  ;  
datafile=['results/2fluoTest/' zname]

load(datafile)

f1=100;
f2=100;
M=f2/f1;

d0=25+zPlanesmm(z);
dLens=f1-d0;

zf=f1;
f=f1;


%estimate detector coords, then referenciate to the center of the detector
detPos=voxelSize.*detMatrix(:,1:3);
maxX=max(detPos(:,1));
maxY=max(detPos(:,2));
detPos(:,1)=detPos(:,1)-maxX/2;
detPos(:,2)=detPos(:,2)-maxY/2;
r=sqrt((maxX/2)^2+(maxY/2)^2);

fprintf('Lens dimensions: %f x %f [mm]\n', maxX, maxY);
fprintf('\t\t\t\t\t %i x %i [px]\n', detPerLine, detPerLine);
fprintf('Lens radius: \t %f mmm\n\n', r);
fprintf('Lens NA:\t\t\t%f  \n', r/f);
fprintf('Focusing at:\t\t\t%f  \n', zf);

detPhotons=length(detPos(detectorCounts.detid));
% detectorCounts.v=cos(detectorCounts.v);
%   detPhotons=10000000;
thetaDev=sqrt(detPos(:,1).^2+detPos(:,2).^2)./f;
rDev=f-sqrt(detPos(:,1).^2+detPos(:,2).^2);
rDev=rDev-max(max(rDev));
rDev=rDev.*0;


%Create a plane for the detector
%Define our camera detector
pxXd=1000;
pxYd=1000;
detPlane=zeros(pxXd,pxYd);
pxSized=0.02; %in mm/px
xd=0:pxXd-1; %x coords
yd=0:pxYd-1; %y coords
xd=xd.*pxSized-pxXd*pxSized/2;
yd=yd.*pxSized-pxYd*pxSized/2;
[Xd,Yd]=meshgrid(xd,yd); %[in mm units]


focusedV=zeros(detPhotons,3);
detPlanePhotons=zeros(detPhotons,2);
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
for i=1:100:detPhotons
    lensPos=rayProp(detPos(detectorCounts.detid(i)-1,:),detectorCounts.v(i,:),dLens);
    
    focusedV(i,:)=rotateVect(lensPos,detectorCounts.v(i,:),-sqrt(lensPos(1).^2+lensPos(2).^2)./f1);
    lensPos=rayProp(lensPos,focusedV(i,:),f1+f2);
    focusedV(i,:)=rotateVect(lensPos, focusedV(i,:),-sqrt(lensPos(1).^2+lensPos(2).^2)./f2);
    
    detPlanePhotons(i,:)=rayProp(lensPos,focusedV(i,:),f2);
    
    detPlanePhotons(i,:)=round(detPlanePhotons(i,:)./pxSized)+pxXd/2+1;
    if detPlanePhotons(i,:)<pxXd & detPlanePhotons(i,:)>0
        detPlane(detPlanePhotons(i,1),detPlanePhotons(i,2))=detPlane(detPlanePhotons(i,1),detPlanePhotons(i,2))+detectorCounts.weights(i)*detectorCounts.w0(i);
        
    end
    
    %Preview of the results
%     if mod(i,100000)==0
%         fprintf('Processed %i out of %i photons\n',i,detPhotons);
%         figure(1)
%         imagesc(xd,yd,detPlane);
%         colorbar
%         drawnow
%         figure(2)
%         imagesc(xd./M,yd./M,log(detPlane+1));
%         colorbar
%         drawnow
%     end
    
end

spimVolume(:,:,z)=detPlane;

% figure,
% imagesc(xd./M,yd./M,detPlane);
% colorbar
% 
% figure,
% imagesc(xd./M,yd./M,log(detPlane+1));
% colorbar


detPlane=detPlane./max(max(detPlane));
%% Now make a beautiful plot
figure,
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(gcf,'color','w');

subplot(1,2,1)
imagesc(xd./M,yd./M,detPlane);
xlabel('x [mm]','Interpreter','latex');
ylabel('y [mm]','Interpreter','latex');
colorbar
cmocean('haline')
pbaspect([1 1 1])
title('Fluorescence image','Interpreter','latex')

subplot(1,2,2)
imagesc(xd./M,yd./M,log10(detPlane+1));
xlabel('x [mm]','Interpreter','latex');
ylabel('y [mm]','Interpreter','latex');
colorbar
cmocean('haline')
 pbaspect([1 1 1])
title('Fluorescence log image','Interpreter','latex')
set(gcf,'position',[10,10,1500,750])

tit=sprintf('Monte Carlo SPIM Simulator [z=%.1f mm - \\mu_a=%6.2f cm^{-1}, \\mu''_s=%6.2f cm^{-1}]\n',zPlanesmm(z),mua,mus);
suptitle(tit);
drawnow
end


fileID = fopen('SPIM_vol.img','w');
fwrite(fileID,single(spimVolume),'single');
fclose(fileID);

