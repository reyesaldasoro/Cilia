cells_0     = imread('C:\Users\sbbk034\OneDrive - City, University of London\Desktop\thumbnail_10xSLC25A18;tm1a;WT;338-0016.jpg');


cells_1     = (cells_0(:,:,2));
%%
cells_2     = imopen(cells_1>(1.4*255*graythresh(cells_1)),strel("disk",4));

cells_3     = bwlabel(cells_2);
cells_3_p   = regionprops(cells_3,'Area');

cells_4     = bwlabel(ismember(cells_3,find( ([cells_3_p.Area]>130).*([cells_3_p.Area]<3330) )  ));
cells_4_p   = regionprops(cells_4,'Area','Circularity','Eccentricity','EquivDiameter','EulerNumber','Perimeter','Solidity');

cells_5     = bwlabel(ismember(cells_4,find( ([cells_4_p.Circularity]>0.7))));
cells_5_p   = regionprops(cells_5,'Area','Circularity','Eccentricity','EquivDiameter','EulerNumber','Perimeter','Solidity');

cells_output=cells_0;
cells_output(:,:,1)=cells_output(:,:,1).*uint8(cells_5==0);
imagesc(cells_output)

%%
imagesc(cells_1.*uint8(imerode(watershed(imfilter(255-cells_1,fspecial('Gaussian',32,13)))>0,ones(3))))
