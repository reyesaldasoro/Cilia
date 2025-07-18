

baseDir     = 'C:\Users\sbbk034\OneDrive - City, University of London\Acad\Research\SGUL_Cilia\TIFFS_2025_07_11';
dir0        = dir(strcat(baseDir,filesep,'*.tif'));

numFiles    = numel(dir0);




%%
for k=87%:numFiles

    currFile    = strcat(baseDir,filesep,dir0(k).name);
    b           = imfinfo(currFile);

    clear a
    figure
    for k2=1:numel(b)
        tempImage           = double(imread(currFile,k2));
        a(:,:,:,k2)         = uint8(295*(tempImage)/max(tempImage(:)));        
    end
    if numel(b)>1
        montage(a)
    else
        imagesc(a)
    end
    colormap hot
    title(dir0(k).name,'Interpreter','none')
end

%%
imagesc(a(:,:,1,6));colorbar

%%
for k=91%:numFiles

    currFile    = strcat(baseDir,filesep,dir0(k).name);

    
    b           = imfinfo(currFile);

    clear a
    for k2=1:3:numel(b)
        %tempImage           = double(imread(currFile,k2));
        a(:,:,1,1+((k2)-1)/3)         = double(imread(currFile,k2+2));  
        a(:,:,2,1+((k2)-1)/3)         = double(imread(currFile,k2+1));  
        a(:,:,3,1+((k2)-1)/3)         = double(imread(currFile,k2+0));  
    end

end
%%
ll=4;
imagesc((a(:,:,2,ll)>200)+(a(:,:,2,ll)>500));colorbar
%%
clear a2
    q=(sort(a(:)));
    thres_99 = q(round(0.995*numel(q)));
    for k3=1:size(a,4)
        a2(:,:,:,k3) = ((a(:,:,:,k3))/thres_99);
    end
%%
    imagesc(a2(:,:,:,2));

    %%
    step=4;
    D1 = isosurface(squeeze(a(1:step:end,1:step:end,3,end:-1:1)),500);
    D2 = patch(D1);
    G1 = isosurface(squeeze(a(1:step:end,1:step:end,2,end:-1:1)),500);
    G2 = patch(G1);


    D2.FaceColor='b';D2.EdgeColor='none';D2.FaceAlpha=0.75;
    G2.FaceColor='g';G2.EdgeColor='none';

    
    lighting gouraud
camlight left

axis ij
axis tight
rotate3d on

%%


[q1,q2] = bwlabeln(squeeze(a(:,:,2,:))>900);
[q3,q4] = bwlabeln(squeeze(a(:,:,2,:))>300);

[q5,q6] = bwlabeln(ismember(q3,unique(q3(q1>0))));
q7=regionprops(q5,squeeze(a(:,:,2,:)),'area',"MeanIntensity","MaxIntensity","MinIntensity");
imagesc((max(q3>0,[],3))+(max(q5>0,[],3)))
%%

imagesc(q5(:,:,7))

%%
    step=2;
    D1 = isosurface(squeeze(a(1:step:end,1:step:end,3,end:-1:1)),500);
    D2 = patch(D1);
    G11 = isosurface(q5(1:step:end,1:step:end,end:-1:1),0.500);
    G21 = patch(G11);
    G21.FaceColor='g';G21.EdgeColor='none';

       D2.FaceColor='b';D2.EdgeColor='none';D2.FaceAlpha=0.75;
    lighting gouraud
camlight left

axis ij
axis tight
rotate3d on