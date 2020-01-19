function [detw] = computeWeight(detectorCounts,prop)
%COMPUTEWEIGHT - One line description of what the function or script performs (H1 line)
%Additional file header info
%
% Syntax:  [weights] = computeWeight(detectorCounts,prop)
%
% Inputs:
%    DETECTORCOUNTS - Description
%    PROP - Description
%
% Outputs:
%    WEIGHTS - Description
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


medianum=size(prop,1);
if(medianum<=1)
    error('empty property list');
end

if(isstruct(detectorCounts))
    if(~isfield(detectorCounts,'w0'))
        detw=ones(size(detectorCounts.ppath,1),1);
    else
        detw=detectorCounts.w0;
    end
    for i=1:medianum-1
        detw=detw.*exp(-prop(i+1,1)*detectorCounts.ppath(:,i));
    end
else
    detectorCounts=detectorCounts';
    if(nargin<3)
        w0=detectorCounts(:,end);
    end
    detw=w0(:);
    if(size(detectorCounts,2)>=2*medianum+1)
        for i=1:medianum-1
            detw=detw.*exp(-prop(i+1,1)*detectorCounts(:,i+medianum));
        end
    else
        for i=1:medianum-1
            detw=detw.*exp(-prop(i+1,1)*detectorCounts(:,i+2));
        end
    end
end
