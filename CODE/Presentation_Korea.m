
clear all
close all



baseDir     = 'C:\Users\sbbk034\OneDrive - City, University of London\Acad\Research\SGUL_Cilia\TIFFS_2025_07_11';
%baseDir     = 'C:\Users\sbbk034\OneDrive - City, University of London\Acad\Research\SGUL_Cilia\TIFFS_2025_07_30';
dir0        = dir(strcat(baseDir,filesep,'*OVER.tif'));

numFiles    = numel(dir0);

%% Read all cells and extract ratios, lengths
jet2                    = [0 0 0;jet];
CalibrationFactor       = 4.8438;
cp                      = cellpose(Model="nuclei");
%%

k=5;
tic
disp(k)
shortName{k}            = dir0(k).name(26:34);
currFile                = strcat(baseDir,filesep,dir0(k).name);
CiliaVolume             = readCilia(currFile);

%%
hot2                        = hot;
green2                      = circshift(hot,1,2);
blue2                       = hot2(:,[3 2 1]);
%% DAPI
imagesc(sum(CiliaVolume(:,:,3,:),4)/255)
colormap(blue2)
 set(gca,'position',[0 0 1 1 ]);axis off
filename = 'Presentation_Korea_1.png';
print('-dpng','-r100',filename)

%% Green
imagesc(sum(CiliaVolume(:,:,2,:),4)/255)
colormap(green2.^0.5)
 set(gca,'position',[0 0 1 1 ]);axis off
filename = 'Presentation_Korea_2.png';
print('-dpng','-r100',filename)

 %% Red
imagesc(sum(CiliaVolume(:,:,1,:),4)/255)
colormap(hot2.^0.5)
 set(gca,'position',[0 0 1 1 ]);axis off
filename = 'Presentation_Korea_3.png';
print('-dpng','-r100',filename)
%%
[CiliaVolume,magnification,calibrationFactor] = readCilia(currFile);
 Output                  = segmentCilia(CiliaVolume,cp,magnification,calibrationFactor);

%%


imagesc(Output.FinalNuclei_MIP)
filename = 'Presentation_Korea_4.png';
print('-dpng','-r100',filename)


%%
imagesc(bwdist(Output.FinalNuclei_MIP))


filename = 'Presentation_Korea_5.png';
print('-dpng','-r100',filename)
%%
figure
imagesc(Output.NucleiRegions)
colormap([0 0 0 ; rand(64,3)])
 set(gca,'position',[0 0 1 1 ]);axis off

filename = 'Presentation_Korea_6.png';
print('-dpng','-r100',filename)

%%
imagesc(repmat(Output.NucleiRegions==10,[1 1 3]).*sum(CiliaVolume(:,:,:,:),4)/255/16)

 set(gca,'position',[0 0 1 1 ]);axis off

filename = 'Presentation_Korea_9.png';
print('-dpng','-r100',filename)
