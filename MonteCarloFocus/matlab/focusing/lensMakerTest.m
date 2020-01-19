clear
clc
close all
% load results/simulation1/detectorCounts_z0mm.mat
% load results/simulation1/simulationData.mat
load results/detectorCounts_z0mm.mat
load results/simulationData.mat

%Lensmaker
M=1;
f=20;
% d0=lensZ;
% d0=f*(1+1/M);
d0=25;
M=f/(d0-f)
zf=1/(1/f-1/d0);

M=1;


if isinf(zf)
    zf=50;
end


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
% detPhotons=1000000;
thetaDev=sqrt(detPos(:,1).^2+detPos(:,2).^2)./f;
rDev=f-sqrt(detPos(:,1).^2+detPos(:,2).^2);
rDev=rDev-max(max(rDev));


%Create a plane for the detector
%Define our camera detector
pxXd=1000;
pxYd=1000;
detPlane=zeros(pxXd,pxYd);
pxSized=0.1; %in mm/px
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

figure(1)
xlabel('mm');
ylabel('mm');
figure(2)
xlabel('mm');
ylabel('mm');
for i=1:1000:detPhotons
    focusedV(i,:)=rotateVectSph(detPos(detectorCounts.detid(i)-1,:),detectorCounts.v(i,:),-thetaDev(detectorCounts.detid(i)-1));
    
    detPlanePhotons(i,:)=rayProp(detPos(detectorCounts.detid(i)-1,:),focusedV(i,:),zf-rDev(detectorCounts.detid(i)-1));
     
    detPlanePhotons(i,:)=round(detPlanePhotons(i,:)./pxSized)+pxXd/2+1;
     if detPlanePhotons(i,:)<pxXd & detPlanePhotons(i,:)>0
        detPlane(detPlanePhotons(i,1),detPlanePhotons(i,2))=detPlane(detPlanePhotons(i,1),detPlanePhotons(i,2))+detectorCounts.weights(i);     
     
     end
    
         if mod(i,100000)==0
        fprintf('Processed %i out of %i photons\n',i,detPhotons);
        figure(1)
        imagesc(xd,yd,detPlane);
        colorbar
        drawnow
        figure(2)
        imagesc(xd./M,yd./M,log(detPlane+1));
        colorbar
        drawnow
    end
  
end


  figure,
    imagesc(xd./M,yd./M,detPlane);
    colorbar
    
        figure,
    imagesc(xd./M,yd./M,log(detPlane+1));
    colorbar


% tic
% figure,
% for i=1:detPhotons
%     focusedV(i,:)=rotateVect(detPos(detectorCounts.detid(i)-1,:),detectorCounts.v(i,:),thetaDev(detectorCounts.detid(i)-1));
%     detPlanePhotons(i,:)=round(rayProp(detPos(detectorCounts.detid(i)-1,:),focusedV(i,:),zf)/pxSize)+pxX/2;
% %     detPlanePhotons(i,:)
%     if detPlanePhotons(i,:)>0 & detPlanePhotons(i,:)<pxX
% %         disp('hit');
%         detPlane(detPlanePhotons(i,1),detPlanePhotons(i,2))= detPlane(detPlanePhotons(i,1),detPlanePhotons(i,2))+detectorCounts.weights(i);
%     else
% %         disp('out')
%     end
%     if mod(i,100000)==0
%         fprintf('Processed %i out of %i photons\n',i,detPhotons);
%         imagesc(log(detPlane));
%         drawnow;
%     end
% end
% toc
% figure,
% imagesc(log(detPlane));


% 
% rad2deg=150/pi;
% 
% h=0:50;  %mm
% 
% R1=50;   %mm
% theta1=h./R1;
% 
% figure,
% plot(h,theta1.*rad2deg);
% xlabel('height')
% ylabel('thetaOut')

