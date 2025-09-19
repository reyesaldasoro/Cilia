function greenPeaks = segmentGreenPeaks(Green_MIP,prom)

% Detect the cilia peaks based on findpeaks in 1D,
[rows,cols ] = size(Green_MIP);
% First start with filtering the data
cilia_MIP_LPF                   = imfilter(Green_MIP,fspecial('gaussian',5,3),'replicate');
greenPeaks_0                    = zeros(size(cilia_MIP_LPF));
% run over all the rows and detect peaks, use a very high prominence to
% only capture significant peaks
for k=1:size(cilia_MIP_LPF,2)
    [~,locPeaks]                = findpeaks(cilia_MIP_LPF(k,:),'MinPeakProminence',prom,'MinPeakDistance',6);
    greenPeaks_0(k,locPeaks)    = 1;
end
% now extend to the regions around those peaks, use Otsu
thres_background                =  max(cilia_MIP_LPF(:))*graythresh(cilia_MIP_LPF/max(cilia_MIP_LPF(:)));
if (sum(cilia_MIP_LPF(:)>(thres_background))/rows/cols)>0.05
    greenPeaks_1                    = bwlabel(cilia_MIP_LPF>(3*thres_background));
else
    greenPeaks_1                    = bwlabel(cilia_MIP_LPF>(thres_background));
end
% if nothing is detected neighbouring the peaks, just dilate a bit
if (sum(greenPeaks_1(:)))==0
    greenPeaks                      = imdilate(greenPeaks_0,strel('disk',2));
else
    greenPeaks_1P                   = regionprops(greenPeaks_1,cilia_MIP_LPF,'area','MeanIntensity','MaxIntensity');
    %keep regions that touch the peaks
    ciliaToKeep                     = unique(greenPeaks_1(greenPeaks_0>0));
    greenPeaks1                     = ismember(greenPeaks_1,ciliaToKeep(2:end));
    %discard regions that are too large (blobs)
    greenPeaks2                     = ismember(greenPeaks_1,find([greenPeaks_1P.Area]<( mean([greenPeaks_1P.Area])+ 3*std([greenPeaks_1P.Area]))));
    greenPeaks                      = greenPeaks2.*greenPeaks1;
end