% clear all;
% clc
% close all;

%Create a plane for the lens [in mm]
pxX=200;
pxY=200;
pxSize=0.1; %in mm/px
lensZ=500; %in px
lensZ=lensZ*pxSize; %in mm

% Define fluorophore
fluoroPos=[100, 100, 0]; %in px units
fluoroPos=[fluoroPos(1)-pxX/2,fluoroPos(2)-pxY/2,fluoroPos(3)];
fluoroPos=fluoroPos*pxSize; %in  mm

%Lensmaker
M=1;
f=200;
d0=lensZ;
d0=f*(1+1/M);
% d0=50
d0=f;
zf=1/(1/f-1/d0);
zf=195;
if isinf(zf)
    zf=50;
end


% detPlanePhotons=zeros(detPhotons,2);

%
fluoroV=zeros(pxX*pxY,3);
fluoroVr=zeros(pxX*pxY,3);
fluoroP=zeros(pxX*pxY,2);
x=0:pxX-1; %x coords
y=0:pxY-1; %y coords
x=x.*pxSize-pxX*pxSize/2;
y=y.*pxSize-pxY*pxSize/2;
[X,Y]=meshgrid(x,y); %[in mm units]

%Array for detected photons
focusedV=zeros(pxX*pxY,5); %[vx,vy,vz,x,y]
 
thetaDev=sqrt(X.^2+Y.^2)./f;
rDev=f-sqrt(f^2-(X.^2+Y.^2));
rDev=rDev-max(max(rDev));
rDev=rDev*0;
% rDev=-rDev;

figure,
hold on,
for i=1:pxX
    for j=1:pxY
        v=[x(i)-fluoroPos(1),y(j)-fluoroPos(2),d0+rDev(i,j)-fluoroPos(3)];
        fluoroV(pxX*(i-1)+j,:)=v./norm(v);
        fluoroVr(pxX*(i-1)+j,:)=rotateVect([x(i), y(j)],v,-thetaDev(i,j));
         fluoroVr(pxX*(i-1)+j,:)=rotateVect([x(i), y(j)],fluoroVr(pxX*(i-1)+j,:),-thetaDev(i,j));
        [zp] = rayProp([x(i), y(j)], fluoroVr(pxX*(i-1)+j,:),zf-rDev(i,j));
        if ~isnan(zp)
            fluoroP(pxX*(i-1)+j,:)=zp;
        end
        %         quiver3(fluoroPos(1),fluoroPos(2),fluoroPos(3),v(1),v(2),v(3));
        if mod(j-1,40)==0 &&  mod(i-1,40)==0
            plot3([x(i),fluoroPos(1)],[y(j),fluoroPos(2)],[d0+rDev(i,j),fluoroPos(3)]);
            plot3([zp(1),x(i)],[zp(2),y(j)],[zf+d0,d0+rDev(i,j)]);
        end
        
    end
    hold on,
end
grid on;
s=surf(X,Y,ones(pxX,pxY).*(d0+rDev),'FaceAlpha',0.5); %Plot the surface
s.EdgeColor = 'none';
xlabel('X');
ylabel('Y');
zlabel('Z');
view(133,-16)
% set(gca,'CameraPosition',[2 2 2]);

%Define our camera detector
pxXd=100;
pxYd=100;
detPlane=zeros(pxXd,pxYd);
pxSized=0.01; %in mm/px

xd=0:pxXd-1; %x coords
yd=0:pxYd-1; %y coords
xd=xd.*pxSized-pxXd*pxSized/2;
yd=yd.*pxSized-pxYd*pxSized/2;
[Xd,Yd]=meshgrid(xd,yd); %[in mm units]
s=surf(Xd,Yd,ones(pxXd,pxYd).*(d0+zf),'FaceAlpha',0.5); %Plot the surface
s.EdgeColor = 'none';

fluoroP=round(fluoroP./pxSized)+pxXd/2+1;


for i=1:length(fluoroP)
    
    if fluoroP(i,:)<pxXd & fluoroP(i,:)>0
        detPlane(fluoroP(i,1),fluoroP(i,2))=detPlane(fluoroP(i,1),fluoroP(i,2))+1;
        
        
    end
    
end
    
    figure,
    imagesc(xd,yd,detPlane);
    colorbar
    
        figure,
    imagesc(xd,yd,log(detPlane+1));
    colorbar
    
    
