function Output = segmentCilia(CiliaVolume,cp,magnification,calibrationFactor)

%% Separate the channels and calculate dimensiones
DAPI                    = squeeze(CiliaVolume(:,:,3,:));
Green                   = squeeze(CiliaVolume(:,:,2,:));
BasalBody               = squeeze(CiliaVolume(:,:,1,:));
BasalBody_MIP           = max(BasalBody,[],3);
DAPI_MIP                = max(DAPI,[],3);
Green_MIP               = max(Green,[],3);
[rows,cols,numSlices]   = size(DAPI);

%% Segment DAPI based on
% (a) CellPose for the maximum intensity projection and 
% (b) Magnification of 060x
if ~exist('cp','var')
    cp                  = cellpose(Model="nuclei");
end

if ~exist('magnification','var')
    magnification       = 60;
end


%NucleiSegmented_MIP     =  segmentCells2D(cp,max(DAPI,[],3),ImageCellDiameter=60);

switch magnification
    case 20
        NucleiSegmented_MIP     =  segmentCells2D(cp,DAPI_MIP,ImageCellDiameter=40, CellThreshold=0,FlowErrorThreshold=1);
        ciliaSize               = 2;
        greenT_H                = 250;
        greenT_L                = 150;
        distFromNucleus         = 8;

    case 40
        NucleiSegmented_MIP     =  segmentCells2D(cp,DAPI_MIP,ImageCellDiameter=70, CellThreshold=0,FlowErrorThreshold=1);
        ciliaSize               = 20;
        greenT_H                = 250;
        greenT_L                = 150;
        distFromNucleus         = 16;

    case 60
        NucleiSegmented_MIP     =  segmentCells2D(cp,DAPI_MIP,ImageCellDiameter=90, CellThreshold=-2,FlowErrorThreshold=1);
        ciliaSize               = 50;
        greenT_H                = 900;
        greenT_L                = 300;
        distFromNucleus         = 24;

    case 100
        NucleiSegmented_MIP     =  segmentCells2D(cp,DAPI_MIP,ImageCellDiameter=120, CellThreshold=-2,FlowErrorThreshold=1);
        ciliaSize               = 70;
        distFromNucleus         = 30;

end
%NucleiSegmented_MIP_P   = regionprops(NucleiSegmented_MIP,'Area','Centroid','BoundingBox','MajorAxisLength','Circularity','Orientation','Eccentricity','MinorAxisLength','Orientation');

% (b) Intensity threshold for the volume
for k=1:numSlices
    tt(k) = 16*255*graythresh(uint8(DAPI(:,:,k)/16));
end
tt2 = sort(tt);

if numSlices>1
    NucleiSegmented         = (smooth3(DAPI,"gaussian",9)>tt2(3)).*repmat(NucleiSegmented_MIP,[1 1 numSlices]);
else
    NucleiSegmented         = (NucleiSegmented_MIP);
    %NucleiSegmented         = (imfilter(DAPI,fspecial("gaussian",9))>tt2(1)).*(NucleiSegmented_MIP);
end
% smooth results
for k=1:numSlices
    NucleiSegmented(:,:,k) = imfill(NucleiSegmented(:,:,k),'holes');
    NucleiSegmented(:,:,k) = imclose(NucleiSegmented(:,:,k),ones(7));
    NucleiSegmented(:,:,k) = imopen(NucleiSegmented(:,:,k),ones(5));
end

%NucleiSegmented_P           =regionprops(NucleiSegmented,'Area','Centroid','BoundingBox','MajorAxisLength','Circularity');
% cellpose for all levels gives some very strange segmentations
 % for k=1:numSlices
 %     disp(k)
 %     NucleiSegmented3(:,:,k)          =  segmentCells2D(cp,DAPI(:,:,k),ImageCellDiameter=90,CellThreshold=0,FlowErrorThreshold=10);
 % end

%%
% Perform a hysteresis thresholding for the cilia
[q1,q2]                 = bwlabeln(Green>greenT_H);
[q3,q4]                 = bwlabeln(Green>greenT_L);
q5                      = ismember(q3,unique(q3(q1>0)));
[CiliaSegmented,q6]     = bwlabeln(ismember(q3,unique(q3(q1>0))));
%CiliaSegmented_P        = regionprops(CiliaSegmented,squeeze(Green),'area',"MeanIntensity","MaxIntensity","MinIntensity");
CiliaSegmented_MIP      = max(CiliaSegmented,[],3);
%CiliaSegmented_MIP_P    = regionprops(CiliaSegmented_MIP,'Area','Centroid','BoundingBox','MajorAxisLength','Circularity','Orientation','Eccentricity','MinorAxisLength','Orientation');

%% Discard cells that touch boundary
AllCells                = unique(NucleiSegmented);
TouchingCells           = union(unique(NucleiSegmented_MIP([1:3 end-3:end],:)),unique(NucleiSegmented_MIP(:,[1:3 end-3:end])));
NucleiSegmented_NotBorder= ismember(NucleiSegmented_MIP,setdiff(AllCells,TouchingCells));
NucleiSegmented_Border  = ismember(NucleiSegmented_MIP,(TouchingCells(2:end)));

%% Discard cilia that are FAR from DAPI and are SMALL
DistFromDAPI            = bwdist(NucleiSegmented_MIP>0);
DistFromDAPI_notTouching = bwdist(NucleiSegmented_NotBorder);


CiliaSegmented_P2       = regionprops(CiliaSegmented,repmat(DistFromDAPI,[1 1 numSlices]),'area',"MeanIntensity","MaxIntensity","MinIntensity");
CiliaSegmented_P3       = regionprops(CiliaSegmented,repmat(DistFromDAPI_notTouching,[1 1 numSlices]),'area',"MeanIntensity","MaxIntensity","MinIntensity");

CiliaToKeep1             = find(([CiliaSegmented_P3.MinIntensity]<distFromNucleus)&([CiliaSegmented_P3.Area]>ciliaSize));
CiliaToKeep             = ismember(CiliaSegmented_MIP,CiliaToKeep1);

%% Determine the final Nuclei and Cilia
% Cilia can be discarded by size or distance to nuclei
% Nuclei can be discarded by contact with boundary
FinalNuclei_MIP             = NucleiSegmented_MIP.*NucleiSegmented_NotBorder;
FinalCilia_MIP              = CiliaSegmented_MIP.*CiliaToKeep;

FinalNuclei_MIP_P           = regionprops(FinalNuclei_MIP,'Area','Centroid','BoundingBox','MajorAxisLength','Circularity','Orientation','Eccentricity','MinorAxisLength','Orientation');
FinalCilia_MIP_P            = regionprops(FinalCilia_MIP,'Area','Centroid','BoundingBox','MajorAxisLength','Circularity','Orientation','Eccentricity','MinorAxisLength','Orientation');



%% Segment Basal Body
% Only consider those regions that are next to Cilia = distance less than
% 10 pixels
distFromCilia               = bwdist(FinalCilia_MIP>0);
BasalBody_MIP               = max(BasalBody,[],3);
BasalBody_ROI               = bwlabel((FinalCilia_MIP==0).*(distFromCilia<=10));
BasalBody_ROI_P             = regionprops(BasalBody_ROI,BasalBody_MIP,'Area','MaxIntensity','MeanIntensity','MinIntensity');

% Only bright pixels, would 150 be a good threshold or determine based on 
meanROI                     = mean(BasalBody_MIP(BasalBody_ROI>0));
stdROI                      = std (BasalBody_MIP(BasalBody_ROI>0));


BasalBody_0                 = ((BasalBody_ROI>0).*(BasalBody_MIP>(meanROI+3*stdROI)));
BasalBody_1                 = imclose(BasalBody_0,ones(5));
BasalBody_1_L               = bwlabel(BasalBody_1);
BasalBody_1_P               = regionprops(BasalBody_1_L,BasalBody_MIP,'Area','MaxIntensity','MeanIntensity','MinIntensity');
BasalBody_2                 = (ismember(BasalBody_1_L,find([BasalBody_1_P.Area]>2)));
%%
BasalBody_3                 = imdilate(BasalBody_2,ones(6,6));
BasalBody_4                 = bwlabel((FinalCilia_MIP==0).*BasalBody_3);


%%
BasalBody_2_P               = regionprops(BasalBody_4,BasalBody_MIP,'Area','MaxIntensity','MeanIntensity','MinIntensity','Centroid');

% 

% if numel(FinalCilia_MIP_P)>0
%     distFromCilia               = bwdist(CiliaSegmented_MIP>0);
%     for k=1:50
%         qq(k)=max(BasalBody_MIP(distFromCilia==k));
%     end
%     plot(qq/(qq(1)))
% end
% 
% imagesc((CiliaSegmented_MIP>0)+(BasalBody_MIP<150).*(distFromCilia<10));colorbar


%% Prepare output

Output.Input_MIP_RGB(:,:,1) = 1.5*BasalBody_MIP;
Output.Input_MIP_RGB(:,:,2) = Green_MIP;
Output.Input_MIP_RGB(:,:,3) = DAPI_MIP;

Output.FinalNuclei_MIP      = FinalNuclei_MIP;
Output.FinalCilia_MIP       = FinalCilia_MIP;
Output.FinalBasalBody_MIP   = BasalBody_4;

Output.FinalNuclei_MIP_P    = FinalNuclei_MIP_P;
Output.FinalCilia_MIP_P     = FinalCilia_MIP_P;
Output.FinalBasalBody_MIP_P = BasalBody_2_P;
Output.CiliaLength          = [Output.FinalCilia_MIP_P.MajorAxisLength]/ calibrationFactor;
Output.NucleiLength         = [Output.FinalNuclei_MIP_P.MajorAxisLength]/ calibrationFactor;


%% output for display, for visualisation cilia and basal body are dilated

BasalBody_3                 = imdilate(BasalBody_4,ones(5));
FinalCilia_MIP_Dil          = imdilate(FinalCilia_MIP,ones(5));

Output.FinalCombination     = 2*(BasalBody_3==0).*(FinalCilia_MIP_Dil>0) + (BasalBody_3==0).*(FinalCilia_MIP_Dil==0).*(Output.FinalNuclei_MIP>0) +3*(BasalBody_3>0);


Output.FinalCombination_RGB(:,:,1) = Output.FinalCombination  ==3 ;
Output.FinalCombination_RGB(:,:,2) = Output.FinalCombination  ==2 ;
Output.FinalCombination_RGB(:,:,3) = Output.FinalCombination  ==1 ;

%imagesc(Output.FinalCombination_RGB)

%%
Output.TotalNuclei          = sum([Output.FinalNuclei_MIP_P.Area]>0);
Output.TotalCilia           = sum([Output.FinalCilia_MIP_P.Area]>0);
Output.TotalBasal           = sum([Output.FinalBasalBody_MIP_P.Area]>0);

Output.Ratio_C_N            = Output.TotalCilia/Output.TotalNuclei;
Output.Ratio_B_C            = Output.TotalBasal/Output.TotalCilia;

Output.BasalBody_MIP        = BasalBody_MIP;
Output.DAPI_MIP             = DAPI_MIP;        
Output.Green_MIP            = Green_MIP;









% k=9;
% imagesc(DAPI(:,:,k).*(NucleiSegmented3(:,:,k)>0))  ;colorbar  
% 
% %%
% 
% k=1;
% imagesc(NucleiSegmented2(:,:,k))  ;colorbar  
