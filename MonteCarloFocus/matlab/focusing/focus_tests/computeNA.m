function [image] = computeNA(detectorCounts,NA,xSize,ySize)
%COMPUTENA - This function computes the information from the detectors
%to obatin an image given a numerical aperture
%
% Syntax:  [image] = computeNA(detectorCounts,NA,xSize,ySize)
%
% Inputs:
%    DETECTORCOUNTS - Structure with the detected photons form MCX
%    NA - Numerical aperture (can be a vector or a single NA)
%    XSIZE - Output image size in the x direction
%    YSIZE - Output image size in the y direction
%
% Outputs:
%    IMAGE - Output image or images (ySize,xSize,number of NA to compute)
%
% Example:
%
% See also:
%
% $Author: Asier Marcos Vidal $    $Date: 18-Oct-2018$    $Revision: 0.1 $
% Copyright: 
%           BiiG - Biomedical Imaging and Instrumentation Group
%           UC3M - Universidad Carlos III de Madrid
%----------------------------- BEGIN CODE ---------------------------------

fprintf('Computing numerical apertures...')
nNA=length(NA);

%initialize output
image=zeros(ySize,xSize,nNA);

%initialize sensor
cameraSensor=zeros(ySize*xSize,1);


%compute detection angles
detectorCounts.ang=asin(detectorCounts.v(:,3));
detectorCounts.detid=detectorCounts.detid-1;
for i=1:nNA
acceptedCounts=detectorCounts.ang>acos(NA(i)/1.33);
detectorCountsNA.detid=detectorCounts.detid(acceptedCounts);
detectorCountsNA.ppath=detectorCounts.ppath(acceptedCounts);
detectorCountsNA.v=detectorCounts.v(acceptedCounts);
detectorCountsNA.ang=detectorCounts.ang(acceptedCounts);
 detectorCountsNA.weights=detectorCounts.weights(acceptedCounts);
% detectorCountsNA.weights=0;

[a,~,c] = unique(detectorCountsNA.detid);
out = [a, accumarray(c,detectorCountsNA.weights)];
cameraSensor(out(:,1))=out(:,2);
image(:,:,i)=reshape(cameraSensor,[ySize,xSize]);
end
