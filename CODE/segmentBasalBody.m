function basalBodyPeaks = segmentBasalBody(BasalBody_MIP,prom)

% Detect the basal body peaks based on findpeaks in 1D,

% First start with filtering the data
BasalBody_MIP_LPF               = imfilter(BasalBody_MIP,fspecial('gaussian',5,3),'replicate');
% run over all the rows and detect peaks, use a very high prominence to
% only capture significant peaks
basalPeaks_0h                    = zeros(size(BasalBody_MIP_LPF));
for k=1:size(BasalBody_MIP_LPF,2)
    [~,locPeaks]                = findpeaks(BasalBody_MIP_LPF(k,:),'MinPeakProminence',prom,'MinPeakDistance',9);
    basalPeaks_0h(k,locPeaks)    = 1;
end
% find also the peaks in the column dimension to improve detection and
% remove ridges
basalPeaks_0v                    = zeros(size(BasalBody_MIP_LPF));
for k=1:size(BasalBody_MIP_LPF,1)
    [~,locPeaks]                = findpeaks(BasalBody_MIP_LPF(:,k),'MinPeakProminence',prom,'MinPeakDistance',9);
    basalPeaks_0v(locPeaks,k)    = 1;
end

basalPeaks_0                    = imdilate(basalPeaks_0h.*basalPeaks_0v,strel("disk",3));

% now extend to the regions around those peaks, use Otsu  OR Dilate and
% select the brightest ones

basalPeaks_0L                   = bwlabel(basalPeaks_0);
basalPeaks_0P                   = regionprops(basalPeaks_0L,BasalBody_MIP_LPF,'area','MeanIntensity','MaxIntensity');
basalBodyPeaks                 = ismember(basalPeaks_0L,find([basalPeaks_0P.MaxIntensity ]>( mean([basalPeaks_0P.MaxIntensity])+ std([basalPeaks_0P.MaxIntensity]))));

    % thres_background                =  max(BasalBody_MIP_LPF(:))*graythresh(BasalBody_MIP_LPF/max(BasalBody_MIP_LPF(:)));
    % basalPeaks_1                    = bwlabel(BasalBody_MIP_LPF>(2.5*thres_background));
    % basalPeaks_1P                   = regionprops(basalPeaks_1,BasalBody_MIP_LPF,'area','MeanIntensity','MaxIntensity');
    % basalBodyPeaks1                 = ismember(basalPeaks_1,unique(basalPeaks_1(basalPeaks_0>0)));
    % basalBodyPeaks2                 = ismember(basalPeaks_1,find([basalPeaks_1P.Area]<( mean([basalPeaks_1P.Area])+ 3*std([basalPeaks_1P.Area]))));
    % basalBodyPeaks                  = basalBodyPeaks2.*basalBodyPeaks1;

%imagesc(q2.*(1-basalPeaks_2))
%%

% q6=bwlabel(basalPeaks_0);
% [q7]=regionprops(q6,BasalBody_MIP_LPF,'area','MinIntensity','MaxIntensity','MeanIntensity','Centroid');
% imagesc(BasalBody_MIP_LPF.*(1-imdilate(ismember(q6,find([q7.Area]>1)),ones(1))))
% 
% 
% %
% k=265;
%     [q3,locPeaks]=findpeaks(BasalBody_MIP_LPF(k,:),'MinPeakProminence',50,'MinPeakDistance',6);
% 
% plot(1:numel(BasalBody_MIP_LPF(k,:)),BasalBody_MIP_LPF(k,:),'b-',locPeaks,q3,'ro')
% grid on 
% axis tight
