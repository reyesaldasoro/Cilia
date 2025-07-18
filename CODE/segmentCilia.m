function CiliaSegmented = segmentCilia(CiliaVolume)


% Perform a hysteresis thresholding for the cilia
[q1,q2] = bwlabeln(squeeze(CiliaVolume(:,:,2,:))>900);
[q3,q4] = bwlabeln(squeeze(CiliaVolume(:,:,2,:))>300);

[CiliaSegmented,q6] = bwlabeln(ismember(q3,unique(q3(q1>0))));
q7=regionprops(q5,squeeze(CiliaVolume(:,:,2,:)),'area',"MeanIntensity","MaxIntensity","MinIntensity");

% Segment DAPI based on intensity


% Discard cilia that is FAR from DAPI

% Segment Basal Body

