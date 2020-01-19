%close all 
clear cfg cfgs
cfg.nphoton=1e7;
cfg.unitinmm=0.1;
dim=200;
cfg.vol=uint8(ones(dim,dim,2*dim+1));
Tflux=zeros(dim,dim,2*dim+1);

cfg.srctype='gaussian';
cfg.srcpos=[0 0 dim];

cfg.srcdir=[1 0 0 100];
cfg.srcparam1=[0.5 0 0];
cfg.gpuid=1;
% cfg.gpuid='11'; % use two GPUs together
cfg.autopilot=1;
cfg.prop=[0 0 1 1;0.1 0.1 0.9 1.37];
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
%cfg.detpos=[50 50 100 50];
cfg.issaveexit=0;
cfg.isreflect=0;
cfg.issrcfrom0=1;
%cfg.maxdetphoton=1e7;

scanPOs=200;
sPos=linspace(0,199,scanPOs);
for i=1:scanPOs  
    cfg.srcpos=[0 sPos(i) dim];
    [flux]=mcxlab(cfg);
    Tflux=Tflux+flux.data*cfg.tstep;
end
Tflux=Tflux/(max(max(max(Tflux))));
%calculate the flux and partial path lengths for the two configurations
fileID=fopen('illuminationVolume.ill','w');
% fwrite(fileID,Nphotons,'uint32');
fwrite(fileID,Tflux,'single');
fclose(fileID);
 
 figure;

imagesc(log10(abs(squeeze(Tflux(:,50,:)))))
axis equal; colorbar
title('pencil beam at volume center');

xAxis=linspace(0,99,100);%.*cfg.unitinmm;
yAxis=linspace(0,99,100);%.*cfg.unitinmm;
zAxis=linspace(0,99,100);%.*cfg.unitinmm;
% [X,Y,Z]=meshgrid(xAxis,yAxis,zAxis);

figure,

   hs=slice(yAxis,xAxis,zAxis,log10(abs(Tflux)),[50],[50],[0,50]);
   set(hs,'linestyle','none');
axis equal;
xlabel('y (mm)');
ylabel('x (mm)');
zlabel('z (mm)');