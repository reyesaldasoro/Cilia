q=imfinfo('C:\Users\sbbk034\OneDrive - City, University of London\Acad\Research\SGUL_Cilia\videos\Live imaging 04092025.tif');
%%
numFrames = size(q,1);
clear c;
for k0=0:24:550 %numFrames
    currFrame = 1+(k0/24);
   disp(k0)
    for k=2:2:12
        a0(:,:,k/2)=double(imread('C:\Users\sbbk034\OneDrive - City, University of London\Acad\Research\SGUL_Cilia\videos\Live imaging 04092025.tif',k0+k-1));
        b0(:,:,k/2)=double(imread('C:\Users\sbbk034\OneDrive - City, University of London\Acad\Research\SGUL_Cilia\videos\Live imaging 04092025.tif',k0+k));
    end

    a1=max(a0,[],3);
    a2=imfilter(a1,fspecial("gaussian",7,3),'replicate');
    a3=a2-min(a2(:));
    a4=a2>36;

    a6=segmentGreenPeaks(a1,5);

    b=b0(:,:,3);

    k1=0.7;

    c(:,:,1,currFrame)=(b)/max(b(:));
    c(:,:,2,currFrame)=c(:,:,1,currFrame).*(1-a6)+a6*100;
    c(:,:,3,currFrame)=(b)/max(b(:));

    imagesc(c(:,:,:,currFrame))
    drawnow
end

montage(c(:,:,:,1:16))