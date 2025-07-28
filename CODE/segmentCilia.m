function CiliaSegmented = segmentCilia(CiliaVolume)

%% Separate the channels

DAPI                    = squeeze(CiliaVolume(:,:,3,:));
Green                   = squeeze(CiliaVolume(:,:,2,:));


[rows,cols,numSlices]   = size(DAPI);

%%
% Perform a hysteresis thresholding for the cilia
[q1,q2]                 = bwlabeln(Green>900);
[q3,q4]                 = bwlabeln(Green>300);
q5                      = ismember(q3,unique(q3(q1>0)));
[CiliaSegmented,q6]     = bwlabeln(ismember(q3,unique(q3(q1>0))));
q7=regionprops(q5,squeeze(CiliaVolume(:,:,2,:)),'area',"MeanIntensity","MaxIntensity","MinIntensity");
%%

%% Segment DAPI based on CellPose
cp                      = cellpose(Model="nuclei");

%NucleiSegmented(:,:,k)          =  segmentCells2D(cp,max(DAPI,[],3),ImageCellDiameter=60);
for k=1:numSlices
    NucleiSegmented(:,:,k)          =  segmentCells2D(cp,DAPI(:,:,k),ImageCellDiameter=60);
end

% Discard cilia that is FAR from DAPI

% Segment Basal Body

