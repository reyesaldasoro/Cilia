%% Cilia project
% Initial assessment of the positions of the cilia
clear all
close all
clc
%% Define locations and folder to read
warning('off')
baseDirectory    = 'C:\Users\sbbk034\OneDrive - City, University of London\Documents\GitHub\Cilia\Data\mmc8';

% read all the folders with data, use wildcards to discard . and ..
dir0            = dir(strcat(baseDirectory,filesep,'*.jpg'));
numFrames       = size(dir0,1);

currImage       = imread(strcat(strcat(baseDirectory,filesep,dir0(187).name)));

%imagesc(currImage)

[rows,cols,levs] = size(currImage);
thresR          = graythresh(currImage(:,:,1));

rr1              = 0+(1:280);
cc1              = 280+(1:340);
rr2              = 361:640;
cc2              = 641:980;
rr3              = 0+(1:280);
cc3              = 641:980;


imagesc(1*currImage(rr1,cc1,:)+0.5*currImage(rr2,cc2,:)+0.5*currImage(rr3,cc3,:))

%%
filtG           = fspecial("gaussian",9,6);
clear F;
for k=202 % :numFrames
    
    currImage = imread(strcat(strcat(baseDirectory,filesep,dir0(k).name)));

    % remove clock and other elements
    
    thresR                          = 85; %255*graythresh(currImage(:,:,1));
    cilia_0                         = imfilter(currImage(rr1,cc1,3),filtG);
    %
    thresG                          = 65; %255*graythresh(currImage(:,:,1));
    cells_0                         = imfilter(currImage(rr3,cc3,2),filtG);
    cells_1                         = (cells_0>thresG);
    [cells_2,numCells]              = bwlabel(cells_1);
    cells_2P                        = regionprops(cells_2,'Area');
    [q1,q2]                         = max([cells_2P.Area]);
    cells_3                         = cells_2==q2;
    %imagesc(cells_3)
    % remove clock and other elements
    cilia_0([1:68],:,:) = 0; 
    cilia_1                         = (cilia_0>thresR).*((currImage(rr2,cc2,1)<220));
    %imagesc(cilia_1)
    composite                         = currImage(rr2,cc2,:);
    
    cilia_2                         = bwmorph(cilia_1,'majority');
    [cilia_3,numObj(k)]                = bwlabel(cilia_2);
    cilia_3P{k}                        = regionprops(cilia_3,'Area','Centroid','BoundingBox','MajorAxisLength','MinorAxisLength','Orientation');
    hold off
    imagesc(composite);
    hold on
    qq = mesh(-0.1+ double(cilia_1));
    qq.EdgeAlpha=0.1;qq.EdgeColor='c';qq.FaceColor='none';
    qq2 = mesh(-0.1+ 0.5*double(cells_3));
    qq2.EdgeAlpha=0.1;qq2.EdgeColor='g';qq2.FaceColor='none';
    view(-40,60)
    axis ([0 340 0 280 0 1])
    drawnow    
    pause(0.1)
    F(k-121) = getframe(gcf);
end
%% For one frame rotate around to show how it matches the still frame
% pan around elevation and azimuth
clear F;
k4=1;
for k1=90:-3:50
    view(0,k1)
    drawnow    
    pause(0.1)
    F(k4) = getframe(gcf);
    k4=k4+1;
end
for k2=[0:5:60 60:-5:-60 -60:5:0]
    view(k2,k1)
    drawnow    
    pause(0.1)
    F(k4) = getframe(gcf);
    k4=k4+1;
end
for k1=50:3:90
    view(0,k1)
    drawnow    
    pause(0.1)
    F(k4) = getframe(gcf);
    k4=k4+1;
end


%% interpolate to discard cases where the decapitated rejoins

%%



   v = VideoWriter('cilia_mmc8_2', 'MPEG-4');
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

        imwrite(imGif,mapGif,'cilia_mmc8_2.gif', 'DelayTime',0,'LoopCount',inf); %g443800