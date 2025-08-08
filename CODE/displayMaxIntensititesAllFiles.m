
clear all
close all



baseDir     = 'C:\Users\sbbk034\OneDrive - City, University of London\Acad\Research\SGUL_Cilia\TIFFS_2025_07_11';
%baseDir     = 'C:\Users\sbbk034\OneDrive - City, University of London\Acad\Research\SGUL_Cilia\TIFFS_2025_07_30';
dir0        = dir(strcat(baseDir,filesep,'*OVER.tif'));

numFiles    = numel(dir0);


%% read all files and display maximum intensity per channel
for k=1:numFiles
    disp(k)
    shortName{k}            = dir0(k).name(26:34);
    currFile                = strcat(baseDir,filesep,dir0(k).name);
    CiliaVolume             = readCilia(currFile);
    maxIntensities(k,:)     = max(max(max((CiliaVolume))),[],4);
    % meanIntensities(k,:)    = mean(mean(mean((CiliaVolume))),4);
    % medianIntensities(k,:)  = median(median(median((CiliaVolume))),4);
    % minIntensities(k,:)     = min(min(min((CiliaVolume))),[],4);
    %max95Intensities(k,:)     =
    MIP_channels            = (max((((CiliaVolume))),[],4));
    basalbodies(:,:,k)      = MIP_channels(:,:,1);
    DAPI(:,:,k)             = MIP_channels(:,:,3);
    Green(:,:,k)            = MIP_channels(:,:,2);
end
%%
hot2                        = hot;
green2                      = circshift(hot,1,2);
blue2                       = hot2(:,[3 2 1]);
figure(1)
montage(uint8(DAPI(:,:,:)/16))
colormap(blue2)
figure(2)
montage(uint8(Green(:,:,:)/4))
colormap(green2)
figure(3)
montage(uint8(basalbodies(:,:,:)/2))
colormap(hot2)
%%
h0=figure;
h1=gca;
hold off
plot(1:numFiles,maxIntensities(:,1),'r-o',1:numFiles,maxIntensities(:,2),'g-o',1:numFiles,maxIntensities(:,3),'b-o')
grid on

% h0=figure;
% h1=gca;
% hold off
% plot(1:numFiles,minIntensities(:,1),'r-o',1:numFiles,minIntensities(:,2),'g-o',1:numFiles,minIntensities(:,3),'b-o')
% grid on
% 
% h0=figure;
% h1=gca;
% hold off
% plot(1:numFiles,meanIntensities(:,1),'r-o',1:numFiles,meanIntensities(:,2),'g-o',1:numFiles,meanIntensities(:,3),'b-o')
% grid on
% 
% h0=figure;
% h1=gca;
% hold off
% plot(1:numFiles,medianIntensities(:,1),'r-o',1:numFiles,medianIntensities(:,2),'g-o',1:numFiles,medianIntensities(:,3),'b-o')
% grid on
title('Maximum Intensities per Channel')
h1.TickLabelInterpreter = 'none';
h1.XTick                = 1:numFiles;
h1.XTickLabel           = shortName;
h1.XTickLabelRotation   = 90;
h1.Position             = [  0.0674    0.27    0.9112    0.60];
h0.Position             = [  643.4000   97.8000  812.0000  333.6000];