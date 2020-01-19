%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIM -  Selective plane illumination microscopy - Monte-Carlo simulator
%
% This simulator uses MCX to simulate light propagation and collection in a
% lightsheet microscope.
%
% Simulations will be runned using 2 configurations, ilcfg and detcfg
%
% by Asier Marcos Vidal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Housekeeping (if needed)
close all;
clear
%% Medium
% Volume size
xSize=200; %in px
ySize=200; %in px
zSize=500; %in px
voxelSize=0.1; %in mm/px

% Optical properties of the medium
mua=0.0025;
mus=1;
g=0.9;
n=1.33;

%Conversion to mm^-1
mua=mua*0.1;
mus=mus*0.1;


%% Parameters
%Zplanes
nzplanes=1;  % number of zplanes to be scanned
zstep=1;    % in voxel units
zcenter=0;  % Central plane to be scanned.

zPlanes= linspace(-(nzplanes-1)/2,(nzplanes-1)/2,nzplanes)*zstep+zcenter;
 zPlanes= linspace(-1,5,7)
 zPlanes=0;
%  zPlanes=5;
zPlanesmm=zPlanes*voxelSize;
zPlanes=zPlanes+zSize/2;


%% Illumination
% Source parameters
% sourcePos=[100 100 0]; % [Y, X, Z]
sourceInitialPos=[0 0 250];
sourceEndPos=[199 0 250];
sourcePositions=200;
sourceDiam=1; %in voxel units
sourceStep=(sourceEndPos-sourceInitialPos)./sourcePositions;

sourcePos=[linspace(sourceInitialPos(1),sourceInitialPos(1)+(sourcePositions-1)*sourceStep(1),sourcePositions);...
    linspace(sourceInitialPos(2),sourceInitialPos(2)+(sourcePositions-1)*sourceStep(2),sourcePositions); ...
    linspace(sourceInitialPos(3),sourceInitialPos(3)+(sourcePositions-1)*sourceStep(3),sourcePositions)];
% sourcePos=[linspace(sourceInitialPos(1),sourceEndPos(1),sourcePositions);...
%     linspace(sourceInitialPos(2),sourceEndPos(2),sourcePositions);...
%     linspace(sourceInitialPos(3),sourceEndPos(3),sourcePositions)];
volume=zeros(ySize,xSize,zSize);

% %% Fluorescence
% 
% % Fluorophores in the medium:
% nFluorophores=36; %distributed along x axis
% fluoroInitialPos=[10 1];% [Y, X]
% fluoroEndPos=[90 200];% [Y, X]
% fluoroPos=round([linspace(fluoroInitialPos(1),fluoroEndPos(1),nFluorophores);...
%     linspace(fluoroInitialPos(2),fluoroEndPos(2),nFluorophores)]);
%    fluoroPos=[5 5 5 5 5 5 10 10 10 10 10 10 20 20 20 20 20 20 50 50 50 50 50 50 100 100 100 100 100 100 175 175 175 175 175 175;
%         5 10 20 50 100 175 5 10 20 50 100 175 5 10 20 50 100 175 5 10 20 50 100 175 5 10 20 50 100 175 5 10 20 50 100 175];
% fluoroPlane=zeros(ySize,xSize);%[Y,X]
% fluoroPlane(fluoroPos(1,:),fluoroPos(2,:))=1;
% % fluoroPlane(100, 100)=1;
% % %  figure,
% %  imshow(fluoroPlane);
% %  title('Fluorophores distribution');
%% Fluorescence

% Fluorophores in the medium:
nFluorophores=2; %distributed along x axis
fluoroInitialPos=[10 1];% [Y, X]
fluoroEndPos=[90 200];% [Y, X]
fluoroPos=round([linspace(fluoroInitialPos(1),fluoroEndPos(1),nFluorophores);...
    linspace(fluoroInitialPos(2),fluoroEndPos(2),nFluorophores)]);
   fluoroPos=[5 5 5 5 5 5 10 10 10 10 10 10 20 20 20 20 20 20 50 50 50 50 50 50 100 100 100 100 100 100 175 175 175 175 175 175;
        5 10 20 50 100 175 5 10 20 50 100 175 5 10 20 50 100 175 5 10 20 50 100 175 5 10 20 50 100 175 5 10 20 50 100 175];
fluoroPlane=zeros(ySize,xSize);%[Y,X]
fluoroPos=[25 175;25 175];
fluoroPlane(fluoroPos(1,:),fluoroPos(2,:))=1;
% fluoroPlane(100, 100)=1;
%  figure,
%  imshow(fluoroPlane);
%  title('Fluorophores distribution');


%% Detection
% Number of detectors is limited by the memory. Several simulations will
% be runned moving the detectors to cover the whole space

%xSize= number of lines (if detector size=1)
%ySize = number of detectors per line (if detector size=1)

detLim=4000; %Maximum number of detectors
detSize=2;   %in voxel units
detcfg.maxdetphoton=1e7;
detcfg.issaveexit=1;


%Divide space in regions
nDetectors=xSize*ySize/detSize^2; %detectors needed
detPerLine=ySize/detSize;
nLines=xSize/detSize;
linesPerChunk=floor(detLim/detPerLine);
nChunks=floor(nLines/linesPerChunk);
lastChunk=nLines-linesPerChunk*nChunks;
if lastChunk>0
    nChunks=nChunks+1;
end

%Generate matrix of detector positions
[Yd,Xd,Zd]=meshgrid(1:detSize:ySize,1:detSize:xSize,zSize);
detMatrix=[Yd(:), Xd(:), Zd(:), repmat(detSize/2,length(Xd(:)),1)];

fprintf('Total detectors used: %d\n',(nChunks*linesPerChunk+lastChunk)*detPerLine);
fprintf('\t Detector divided in %d regions\n',nChunks);


%% Prepare Monte-Carlo illumination configuration
simulateIllumination=0;
%Photons
nphotons=1e8;
illcfg.nphoton=nphotons;
%Medium
illcfg.vol=uint8(ones(xSize,ySize,zSize));
illcfg.prop=[0 0 1 1;mua mus g n];  %[mua, mus, g, n]
illcfg.isreflect=0;
illcfg.isrefint=0;

% Source Description
illcfg.srctype='gaussian';
% illcfg.srcpos=sourcePos;
illcfg.srcdir=[0 1 0]; %[Y, X, Z]
illcfg.srcparam1=[sourceDiam/2 0 0 0]; %Size in voxels
illcfg.srcparam2=[0 0 0 0];

% GPU
illcfg.gpuid=1;
illcfg.autopilot=1;
illcfg.isgpuinfo=0;

%Timing
illcfg.tstart=0;
illcfg.tend=5e-9;
illcfg.tstep=5e-9;
illcfg.seed=99999;
illcfg.unitinmm=voxelSize;

%Output type
illcfg.outputtype='flux';

%% Initialize figure
figure(1),
subplot(1,2,1)

xAxis=linspace(-(xSize*voxelSize/2-voxelSize),xSize*voxelSize/2,xSize);
yAxis=linspace(-(ySize*voxelSize/2-voxelSize),ySize*voxelSize/2,ySize);
zAxis=linspace(-(zSize*voxelSize/2-voxelSize),zSize*voxelSize/2,zSize);
% xAxis=linspace(-xSize*voxelSize/2,xSize*voxelSize/2-voxelSize,xSize);
% yAxis=linspace(-ySize*voxelSize/2,ySize*voxelSize/2-voxelSize,ySize);
% zAxis=linspace(0,zSize*voxelSize-voxelSize,zSize);
[Y,X,Z]=meshgrid(yAxis,xAxis,zAxis);
% subplot(223);
%% Run illumination simulation
if simulateIllumination
illvolume=zeros(ySize,xSize,zSize);

% hs=slice(Y,X,Z,log10(abs(double(fcw))),[yAxis(sourcePos(1))],[xAxis(sourcePos(2))],[zAxis(sourcePos(3)+1)]);
hs=slice(Y,X,Z,illvolume,[0],[0],[0]);
set(hs,'linestyle','none');
axis equal;
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
title(sprintf('Simulating illumination... Source %d of %d',1,sourcePositions));
hold on
drawnow
for i=1:sourcePositions
    fprintf('Simulating source position %i of %i\n',i,sourcePositions);
    illcfg.srcpos=sourcePos(:,i);
    flux=mcxlab(illcfg);
    illvolume=flux.data*illcfg.tstep+illvolume;
    hs=slice(Y,X,Z,log10(abs(double(illvolume))),[0],[0],[0]);
    set(hs,'linestyle','none');
    title(sprintf('Simulating illumination... Source %d of %d',i,sourcePositions));
    colorbar
    drawnow
end

%Now normalize output volume to obtain fluorescence excitation

% % Normalize illumination volume
 illvolume=double(illvolume)./max(max(max(illvolume)));
else
   load illvolume
   hs=slice(Y,X,Z,log10(abs(double(illvolume))),[0],[0],[0]);
   set(hs,'linestyle','none');
axis equal;
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
end
title('Illumination')

%% Prepare Monte-Carlo detection configuration

detcfg.nphoton=nphotons;
%Medium
detcfg.vol=uint8(ones(xSize,ySize,zSize));
detcfg.prop=[0 0 1 1;mua mus g n];  %[mua, mus, g, n]
detcfg.isreflect=0;
detcfg.isrefint=0;
detcfg.issrcfrom0=1;

% Source Description
 detcfg.srctype='isotropic';
% illcfg.srcpos=sourcePos;
detcfg.srcdir=[0 1 0]; %[Y, X, Z]
detcfg.srcparam1=[0.5 0 0 0]; %Size in voxels
detcfg.srcparam2=[0 0 0 0];

% GPU
detcfg.gpuid=1;
detcfg.autopilot=1;
detcfg.isgpuinfo=0;

%Timing
detcfg.tstart=0;
detcfg.tend=5e-9;
detcfg.tstep=5e-9;
detcfg.seed=99999;
detcfg.unitinmm=voxelSize;

%Output type
detcfg.outputtype='flux';


%Here we save the configuration

 save('results\\simulationData.mat')

%% Run detection simulation

%Hay que normalizar la energía de excitación
for j=1:nzplanes
detvolume=zeros(ySize,xSize,zSize);
figure(1),
subplot(1,2,2)
hs=slice(Y,X,Z,detvolume,[0],[0],[0]);
set(hs,'linestyle','none');
axis equal;
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');

hold on

maxDet=0;
detectorCounts.detid=[];
detectorCounts.v=[];
detectorCounts.ppath=[];
detectorCounts.w0=[];
scatter3(yAxis(Yd(:)),xAxis(Xd(:)),zAxis(Zd(:)),'r');
for f=1:nFluorophores
    title(sprintf('Simulating detection... Fluorophore %d of %d',f,nFluorophores));
    fprintf('Fluorophore %i of %i\n',i,nChunks);
    detcfg.srcpos=[fluoroPos(1,f),fluoroPos(2,f) zSize/2 ];
    maxDet=0;
    fluoExcitation=illvolume(detcfg.srcpos(1),detcfg.srcpos(2),zPlanes(j));
    for i=1:nChunks
        fprintf('Simulating detector sector %i of %i\n',i,nChunks);
        minDet=maxDet+1;
        maxDet=i*detPerLine*linesPerChunk;
        if maxDet>nDetectors
            maxDet=nDetectors;
        end
        
        detcfg.detpos=detMatrix(minDet:maxDet,:);
        hold on
        if exist('d','var')
            delete(d);
        end
        d=scatter3(yAxis(Yd(1:maxDet)),xAxis(Xd(1:maxDet)),zAxis(Zd(1:maxDet)),'g');
        drawnow
        [flux, detp1]=mcxlab(detcfg);
        if isfield(detp1,'detid')
            detectorCounts.detid=[detectorCounts.detid; detp1.detid+minDet];
            detectorCounts.v=[detectorCounts.v; detp1.v];
            detectorCounts.ppath=[detectorCounts.ppath; detp1.p];
            detectorCounts.w0=[detectorCounts.w0;repmat(fluoExcitation,length(detp1.detid),1)];
        end
                 
    end
%     fcw=fcw.*detvolume(cfg.srcpos(1),cfg.srcpos(2),cfg.srcpos(3));


    detvolume=flux.data*detcfg.tstep.*fluoExcitation+detvolume;
    hs=slice(Y,X,Z,log10(abs(double(detvolume))),[0],[0],[0]);
    set(hs,'linestyle','none');
    colorbar
end

title('Detection');

% clear d hs dp1 f i detLim 


%% Data processing
detectorCounts.weights = computeWeight(detectorCounts,detcfg.prop);

%Here we save the counts with the scanned plane
planeName=sprintf('results\\detectorCounts_z%.1gmm.mat',zPlanesmm(j));
save(char(planeName),'detectorCounts','detvolume','-v7.3');

% 

% NA=[0.2, 0.04, 0.025, 0.01, 0.005, 0.001];
% 
% image = computeNA(detectorCounts,NA,xSize/detSize,ySize/detSize);
% figure,
% nNA=length(NA);
% for i=1:nNA
%     subplot(nNA/2,nNA/2,i)
%     imagesc([xAxis(1),xAxis(end)],[yAxis(1),yAxis(end)],image(:,:,i)')
%     axis equal
%     xlabel('x (mm)');
%     ylabel('y (mm)');
%     title(sprintf('NA %.2g',NA(i)))   
% end
end



% 
% 
% 
% el=asin(detectorCounts.v(:,3));  % elevation angle of v, el=0 parallel to surface, el=pi/2 at normal dir
% 
% edges=linspace(min(el),max(el),100);  % angle bins for hisotogram
% 
% [ct, bin]=histc(el,edges); % count of photons per angle bin
% 
% % R=cfg.detpos(1,4);         % radius of the det
% R=1;
% hedges=abs(R*sin(edges));  % height of each spherical segment for each bin
% zonearea=2*pi*R*(diff(hedges)); % area of each spherical segment
% 
% detweight=mmcdetweight(detectorCounts, detcfg.prop); % get detected photon weight
% angularweight=accumarray(bin,detweight); % sum total weight per zone
% 
% angularflux=angularweight(1:end-1)'./zonearea;  % calculate flux per angle
% figure,
% polar((edges(1:end-1)+edges(2:end))*0.5, angularflux); % plot flux vs angle
% hold on;
% polar(pi-(edges(1:end-1)+edges(2:end))*0.5, angularflux); % mirror to form a circle

%
%
% %% Prepare Monte-Carlo Simulation
% clear cfg;
%
% %Photons
% nphotons=1e8;
% cfg.nphoton=nphotons;
% %Medium
% cfg.vol=uint8(ones(xSize,ySize,zSize));
% cfg.prop=[0 0 1 1;mua mus g n];  %[mua, mus, g, n]
% cfg.isreflect=0;
% cfg.isrefint=0;
%
% % Source Description
% cfg.srctype='gaussian';
% cfg.srcpos=sourcePos;
% cfg.srcdir=[0 1 0]; %[Y, X, Z]
% cfg.srcparam1=[0.5 0 0 0]; %Size in voxels
% cfg.srcparam2=[0 0 0 0];
%
% % GPU
% cfg.gpuid=1;
% cfg.autopilot=1;
% cfg.isgpuinfo=0;
%
% %Timing
% cfg.tstart=0;
% cfg.tend=5e-9;
% cfg.tstep=5e-9;
% cfg.seed=99999;
% cfg.unitinmm=voxelSize;
%
% %Output type
% cfg.outputtype='energy';
%
% %All parameters stay the same except from source position
%
% for i=1:sourcePositions
%     fprintf('Simulating source position %i of %i\n',i,sourcePositions);
%     cfg.srcpos=sourcePos(:,i);
%     flux=mcxlab(cfg);
%     fcw=flux.data;
%     fcw=flux.data*cfg.tstep; %Here convert from fluence rate to fluence
%     volume=fcw+volume;
% end
% fcw=volume;
% %Prepare Plot
% xAxis=linspace(-(xSize*voxelSize/2-voxelSize),xSize*voxelSize/2,xSize);
% yAxis=linspace(-(ySize*voxelSize/2-voxelSize),ySize*voxelSize/2,ySize);
% zAxis=linspace(-(zSize*voxelSize/2-voxelSize),zSize*voxelSize/2,zSize);
% % xAxis=linspace(-xSize*voxelSize/2,xSize*voxelSize/2-voxelSize,xSize);
% % yAxis=linspace(-ySize*voxelSize/2,ySize*voxelSize/2-voxelSize,ySize);
% % zAxis=linspace(0,zSize*voxelSize-voxelSize,zSize);
% [Y,X,Z]=meshgrid(yAxis,xAxis,zAxis);
% % subplot(223);
% figure,
% % hs=slice(Y,X,Z,log10(abs(double(fcw))),[yAxis(sourcePos(1))],[xAxis(sourcePos(2))],[zAxis(sourcePos(3)+1)]);
% hs=slice(Y,X,Z,log10(abs(double(fcw))),[0],[0],[0]);
%
% % hs=slice(Y,X,Z,abs(double(fcw)),[],[yAxis(sourcePos(2)+1)],[zAxis(sourcePos(3)+1)]);
%
% % hs=slice(log10(abs(double(fcw))),[],[10 45],[1]);
% set(hs,'linestyle','none');
% axis equal; colorbar
% xlabel('x (mm)');
% ylabel('y (mm)');
% zlabel('z (mm)');
% title('a slit source');
%
% %% Simulate fluorescence
%
% % Normalize illumination volume
% vokume=double(volume)./max(max(max(volume)));
%
% % Configure Monte-Carlo Simulation
% % Source Description
% cfg.srctype='isotropic';
% % cfg.srcpos=sourcePos;
% cfg.srcdir=[0 0 1]; %[Y, X, Z]
% cfg.srcparam1=[1 0 0 0]; %Size in voxels
% cfg.srcparam2=[0 0 0 0];
%
% fluoroVolume=zeros(ySize,xSize,zSize);
%
% %By now only one plane
% zSource=sourceInitialPos(3);
% for i=1:nFluorophores
%     %     fprintf('Simulating source position %i of %i\n',i,sourcePositions);
%
%     % Configure source position
%     cfg.srcpos=[fluoroPos(:,i);zSource];
%     % get fluorescence at that position
%     % fcw(cfg.srcpos)
%     % cfg.nphoton=fcw(cfg.srcpos)*nphotons;
%     flux=mcxlab(cfg);
%     fcw=flux.data./max(max(max(flux.data)));
%     fcw=fcw.*volume(cfg.srcpos(1),cfg.srcpos(2),cfg.srcpos(3)); %Here convert from fluence rate to fluence
%     fcw(cfg.srcpos(1),cfg.srcpos(2),cfg.srcpos(3))
%     fluoroVolume=fcw+fluoroVolume;
% end
% figure,
% % hs=slice(Y,X,Z,log10(abs(double(fcw))),[yAxis(sourcePos(1))],[xAxis(sourcePos(2))],[zAxis(sourcePos(3)+1)]);
% hs=slice(Y,X,Z,log10(abs(double(fluoroVolume))),[0],[0],[0,5]);
%
% % hs=slice(Y,X,Z,abs(double(fcw)),[],[yAxis(sourcePos(2)+1)],[zAxis(sourcePos(3)+1)]);
%
% % hs=slice(log10(abs(double(fcw))),[],[10 45],[1]);
% set(hs,'linestyle','none');
% axis equal; colorbar
% xlabel('x (mm)');


% ylabel('y (mm)');
% zlabel('z (mm)');
% title('a slit source');
