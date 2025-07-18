function CiliaVolume = readCilia(currFile)

b                                       = imfinfo(currFile);
numSlices                               = numel(b);

CiliaVolume(b(1).Width,b(1).Height,3,numSlices/3) = 0;
for k2=1:3:numel(b)
    %tempImage           = double(imread(currFile,k2));
    CiliaVolume(:,:,1,1+((k2)-1)/3)     = double(imread(currFile,k2+2));
    CiliaVolume(:,:,2,1+((k2)-1)/3)     = double(imread(currFile,k2+1));
    CiliaVolume(:,:,3,1+((k2)-1)/3)     = double(imread(currFile,k2+0));
end