%% Cilia project
% Initial assessment of the positions of the cilia
clear all
close all
clc
%% Define locations and folder to read
warning('off')
baseDirectory    = 'C:\Users\sbbk034\OneDrive - City, University of London\Documents\GitHub\Cilia\Data\mmc4';

% read all the folders with data, use wildcards to discard . and ..
dir0            = dir(strcat(baseDirectory,filesep,'*.jpg'));
numFrames       = size(dir0,1);

currImage       = imread(strcat(strcat(baseDirectory,filesep,dir0(1).name)));
[rows,cols,levs] = size(currImage);
thresR          = graythresh(currImage(:,:,1));
%%
filtG           = fspecial("gaussian",9,6);

for k=1:numFrames
    currImage = imread(strcat(strcat(baseDirectory,filesep,dir0(k).name)));

    % remove clock and other elements
    currImage([1:100 350:rows],:,:) = 0; 
    currImage(:,[1:100 650:cols],:) = 0; 
    thresR                          = 255*graythresh(currImage(:,:,1));
    cilia_0                         = imfilter(currImage(:,:,1),filtG);
    cilia_1                         = cilia_0>thresR;
    cilia_2                         = bwmorph(cilia_1,'majority');
    [cilia_3,numObj(k)]                = bwlabel(cilia_2);
    cilia_3P{k}                        = regionprops(cilia_3,'Area','Centroid','BoundingBox','MajorAxisLength','Orientation');
    %imagesc(cilia_3);colorbar
    %drawnow
end

%% interpolate to discard cases where the decapitated rejoins
centroids_2=[];
centroids_1=[];
for k=1:numFrames
    centroids_1(k,:) = [cilia_3P{k}(1).Centroid cilia_3P{k}(1).Orientation cilia_3P{k}(1).MajorAxisLength] ;
    if numObj(k)>1
        centroids_2(k,:) = [cilia_3P{k}(2).Centroid  cilia_3P{k}(1).Orientation cilia_3P{k}(1).MajorAxisLength];
    end
end
%%
for k=2:numFrames-1
    centroids_1(k,:) = median([centroids_1(k-1:k+1,:)]);
    centroids_2(k,:) = median([centroids_2(k-1:k+1,:)]);
end
%%

for k=8%:numFrames
    currImage = imread(strcat(strcat(baseDirectory,filesep,dir0(k).name)));

    % remove clock and other elements
    currImage([1:100 350:rows],:,:) = 0; 
    currImage(:,[1:100 650:cols],:) = 0; 
    thresR                          = 255*graythresh(currImage(:,:,1));
    cilia_0                         = imfilter(currImage(:,:,1),filtG);
    cilia_1                         = cilia_0>thresR;
    cilia_2                         = bwmorph(cilia_1,'majority');
    [cilia_3,numObj(k)]                = bwlabel(cilia_2);
    hold on
    imagesc(cilia_3);colorbar

    plot(centroids_1(k,1),centroids_1(k,2),'ro')
    plot(centroids_2(k,1),centroids_2(k,2),'md')
    drawnow
    hold off
end


