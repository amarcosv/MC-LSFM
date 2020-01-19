clear
load simdata
cameraSensor=zeros(nDetectors,1);

% cameraSensor=detectorCounts.detid

[a,~,c] = unique(detectorCounts.detid);
out = [a, accumarray(c,detectorCounts.ppath)];
out(:,1)=out(:,1)-1;
cameraSensor(out(:,1))=double(out(:,2));

cameraImage=reshape(cameraSensor,[ySize/detSize,xSize/detSize]);
figure,

imagesc(cameraImage)
colorbar