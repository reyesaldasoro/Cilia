%% Cilia project
% Initial assessment of the positions of the cilia
clear all
close all
clc
%% Define locations and folder to read
warning('off')
if strcmp('/',filesep)
    baseDirectory    = '/Users/ccr22/Academic/GitHub/Cilia/Data/mmc8';    
else
    baseDirectory    = 'C:\Users\sbbk034\OneDrive - City, University of London\Documents\GitHub\Cilia\Data\mmc8';
end
% read all the folders with data, use wildcards to discard . and ..
dir0                        = dir(strcat(baseDirectory,filesep,'*.jpg'));
numFrames                   = size(dir0,1);

currImage                   = imread(strcat(strcat(baseDirectory,filesep,dir0(1).name)));
[rows,cols,levs]            = size(currImage);
thresR                      = graythresh(currImage(:,:,1));
%%
rr_panel_11                 = 68:360;
cc_panel_11                 = 80:642;



filtG                       = fspecial("gaussian",9,6);

for k=122:numFrames
    currImage               = imread(strcat(strcat(baseDirectory,filesep,dir0(k).name)));
    currImage1              = currImage(rr_panel_11,cc_panel_11,3);
    thresR                  = 255*graythresh(currImage1);
    %currImage2              = currImage1.*uint8(currImage1>thresR);
    %currImage2(currImage2==0)=thresR;
    currImage3              = bwlabel(currImage1>80);
    currImage4              = bwmorph(currImage1>135,'majority');
    keepArea                = (unique(currImage3.*currImage4));
    
    currImage5              = ismember(currImage3,keepArea(2:end));
    
    imagesc(currImage5);colorbar
    
    %caxis([0 255])
%     % remove clock and other elements
%     currImage([1:100 350:rows],:,:) = 0; 
%     currImage(:,[1:100 650:cols],:) = 0; 
%     thresR                          = 255*graythresh(currImage(:,:,1));
%     cilia_0                         = imfilter(currImage(:,:,1),filtG);
%     cilia_1                         = cilia_0>thresR;
%     cilia_2                         = bwmorph(cilia_1,'majority');
%     [cilia_3,numObj(k)]                = bwlabel(cilia_2);
%     cilia_3P{k}                        = regionprops(cilia_3,'Area','Centroid','BoundingBox','MajorAxisLength','MinorAxisLength','Orientation');
    %imagesc(cilia_3)
    pause(0.1)
    drawnow
end

%% interpolate to discard cases where the decapitated rejoins
centroids_2=[];
centroids_1=[];
for k=1:numFrames
    centroids_1(k,:) = [cilia_3P{k}(1).Centroid cilia_3P{k}(1).Orientation cilia_3P{k}(1).MajorAxisLength cilia_3P{k}(1).MinorAxisLength] ;
    if numObj(k)>1
        centroids_2(k,:) = [cilia_3P{k}(2).Centroid  cilia_3P{k}(1).Orientation cilia_3P{k}(2).MajorAxisLength  cilia_3P{k}(2).MinorAxisLength];
    end
end
%%
for k=2:numFrames-1
    centroids_1(k,:) = median([centroids_1(k-1:k+1,:)]);
    centroids_2(k,:) = median([centroids_2(k-1:k+1,:)]);
end
%%
clear F;
for k=2:numFrames
    currImage = imread(strcat(strcat(baseDirectory,filesep,dir0(k).name)));

    % remove clock and other elements
    currImage([1:100 350:rows],:,:) = 0; 
    currImage(:,[1:100 650:cols],:) = 0; 
    thresR                          = 255*graythresh(currImage(:,:,1));
    cilia_0                         = imfilter(currImage(:,:,1),filtG);
    cilia_1                         = cilia_0>thresR;
    cilia_2                         = bwmorph(cilia_1,'majority');
    [cilia_3,numObj(k)]                = bwlabel(cilia_2);

    imagesc(currImage);%colorbar
    hold on
    text(centroids_1(k,1),centroids_1(k,2)+25,2,num2str(centroids_1(k,4),3),'fontsize',22,'color','m')
    text(centroids_2(k,1),centroids_2(k,2)+10,2.3,num2str(centroids_2(k,4),3),'fontsize',22,'color','b')
    %plot(centroids_2(k,1),centroids_2(k,2),'md')
    [X1,Y1,Z1] = ellipsoid(centroids_1(k,1),centroids_1(k,2), 1,centroids_1(k,5)/2,centroids_1(k,4)/2,0.1);
    [X2,Y2,Z2] = ellipsoid(centroids_2(k,1),centroids_2(k,2), 1,centroids_2(k,5)/2,centroids_2(k,4)/2,0.1);
    e1 = surf(X1,Y1,Z1);
    rotate(e1,[0 0 1],(sign(centroids_1(k,3)))*(90-abs(centroids_1(k,3))))
    e2 = surf(X2,Y2,Z2);
    e1.EdgeColor='m';
    e1.EdgeAlpha=0.25;
    e2.EdgeAlpha=0.25;
    e2.EdgeColor='c';
    e1.FaceColor='none';
    e2.FaceColor='none';
    view(50,70)
   axis([350 650 50 350])
    drawnow
    pause(0.1)
    F(k-1) = getframe(gcf);
    hold off
end


%%
k2=k;
for k=50:2:87

     view(k,70)
    drawnow
    pause(0.1)
    F(k2) = getframe(gcf);
    k2=k2+1;
end 
%%
for k=[70:2:90 90:-2:70 ]%30 30:2:90 90:-2:60]
    view(87,k)
    drawnow
    pause(0.1)
    F(k2) = getframe(gcf);
        k2=k2+1;
end 

%%



   v = VideoWriter('cilia_2A', 'MPEG-4');
            open(v);
            writeVideo(v,F);
            close(v);

%% save the movie as a GIF
    [imGif,mapGif] = rgb2ind(F(1).cdata,256,'nodither');
    numFrames = size(F,2);

    imGif(1,1,1,numFrames) = 0;
    for k = 2:numFrames 
      imGif(:,:,1,k) = rgb2ind(F(k).cdata,mapGif,'nodither');
    end

        imwrite(imGif,mapGif,'cilia_2B.gif', 'DelayTime',0,'LoopCount',inf); %g443800