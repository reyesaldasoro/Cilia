function ciliaPeaks = segmentCiliaPeaks(Green_MIP,prom)

% Detect the cilia peaks based on findpeaks in 1D,

% First start with filtering the data
cilia_MIP_LPF                   = imfilter(Green_MIP,fspecial('gaussian',5,3),'replicate');
ciliaPeaks_0                    = zeros(size(cilia_MIP_LPF));
% run over all the rows and detect peaks, use a very high prominence to
% only capture significant peaks
for k=1:size(cilia_MIP_LPF,2)
    [~,locPeaks]                = findpeaks(cilia_MIP_LPF(k,:),'MinPeakProminence',prom,'MinPeakDistance',6);
    ciliaPeaks_0(k,locPeaks)    = 1;
end
% now extend to the regions around those peaks, use Otsu
thres_background                =  max(cilia_MIP_LPF(:))*graythresh(cilia_MIP_LPF/max(cilia_MIP_LPF(:)));
ciliaPeaks_1                    = bwlabel(cilia_MIP_LPF>(thres_background));
ciliaPeaks_1P                   = regionprops(ciliaPeaks_1,cilia_MIP_LPF,'area','MeanIntensity','MaxIntensity');
%keep regions that touch the peaks
ciliaToKeep                     = unique(ciliaPeaks_1(ciliaPeaks_0>0));
ciliaPeaks1                     = ismember(ciliaPeaks_1,ciliaToKeep(2:end));
%discard regions that are too large (blobs)
ciliaPeaks2                     = ismember(ciliaPeaks_1,find([ciliaPeaks_1P.Area]<( mean([ciliaPeaks_1P.Area])+ 3*std([ciliaPeaks_1P.Area]))));
ciliaPeaks                      = ciliaPeaks2.*ciliaPeaks1;