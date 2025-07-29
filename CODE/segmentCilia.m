function CiliaSegmented = segmentCilia(CiliaVolume)

%% Separate the channels and calculate dimensiones
DAPI                    = squeeze(CiliaVolume(:,:,3,:));
Green                   = squeeze(CiliaVolume(:,:,2,:));
[rows,cols,numSlices]   = size(DAPI);

%% Segment DAPI based on
% (a) CellPose for the maximum intensity projection and 

if ~exist('cp','var')
    cp                  = cellpose(Model="nuclei");
end
NucleiSegmented_MIP     =  segmentCells2D(cp,max(DAPI,[],3),ImageCellDiameter=60);
NucleiSegmented_MIP_P   = regionprops(NucleiSegmented_MIP,'Area','Centroid','BoundingBox','MajorAxisLength','Circularity','Orientation','Eccentricity','MinorAxisLength','Orientation');

% (b) Intensity threshold for the volume
for k=1:numSlices
    tt(k) = 16*255*graythresh(uint8(DAPI(:,:,k)/16));
end
tt2 = sort(tt);
NucleiSegmented         = (smooth3(DAPI,"gaussian",9)>tt2(3)).*repmat(NucleiSegmented_MIP,[1 1 numSlices]);

% smooth results
for k=1:numSlices
    NucleiSegmented(:,:,k) = imfill(NucleiSegmented(:,:,k),'holes');
    NucleiSegmented(:,:,k) = imclose(NucleiSegmented(:,:,k),ones(7));
    NucleiSegmented(:,:,k) = imopen(NucleiSegmented(:,:,k),ones(5));
end

NucleiSegmented_P           =regionprops(NucleiSegmented,'Area','Centroid','BoundingBox','MajorAxisLength','Circularity');
% cellpose for all levels gives some very strange segmentations
 % for k=1:numSlices
 %     disp(k)
 %     NucleiSegmented3(:,:,k)          =  segmentCells2D(cp,DAPI(:,:,k),ImageCellDiameter=90,CellThreshold=0,FlowErrorThreshold=10);
 % end

%%
% Perform a hysteresis thresholding for the cilia
[q1,q2]                 = bwlabeln(Green>900);
[q3,q4]                 = bwlabeln(Green>300);
q5                      = ismember(q3,unique(q3(q1>0)));
[CiliaSegmented,q6]     = bwlabeln(ismember(q3,unique(q3(q1>0))));
CiliaSegmented_P        = regionprops(CiliaSegmented,squeeze(Green),'area',"MeanIntensity","MaxIntensity","MinIntensity");
CiliaSegmented_MIP_P    = regionprops(max(CiliaSegmented,[],3),'Area','Centroid','BoundingBox','MajorAxisLength','Circularity','Orientation','Eccentricity','MinorAxisLength','Orientation');

%% Discard cells that touch boundary
AllCells                = unique(NucleiSegmented);
TouchingCells           = union(unique(NucleiSegmented_MIP([1:3 end-3:end],:)),unique(NucleiSegmented_MIP(:,[1:3 end-3:end])));
NucleiSegmented_NotBorder=ismember(NucleiSegmented_MIP,setdiff(AllCells,TouchingCells));

%%
% Discard cilia that is FAR from DAPI
DistFromDAPI            = bwdist(NucleiSegmented_MIP>0);
DistFromDAPI_notToching = bwdist(NucleiSegmented_NotBorder);


CiliaSegmented_P2       = regionprops(CiliaSegmented,repmat(DistFromDAPI,[1 1 numSlices]),'area',"MeanIntensity","MaxIntensity","MinIntensity");
CiliaSegmented_P3       = regionprops(CiliaSegmented,repmat(DistFromDAPI_notToching,[1 1 numSlices]),'area',"MeanIntensity","MaxIntensity","MinIntensity");

%% Segment Basal Body
% k=9;
% imagesc(DAPI(:,:,k).*(NucleiSegmented3(:,:,k)>0))  ;colorbar  
% 
% %%
% 
% k=1;
% imagesc(NucleiSegmented2(:,:,k))  ;colorbar  
