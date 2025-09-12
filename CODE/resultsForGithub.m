hot2=hot;
blue3=hot2(:,[3 2 1]);
green3=hot2(:,[2 1 3 ]);
colormap(green3)


imagesc(Output.Input_MIP_RGB/255)
filename = '../Figures/Figure_1_RGB_MIP.png';

imagesc(max(CiliaVolume(:,:,1,:),[],4));
colormap(hot.^0.5)
filename = '../Figures/Figure_1_Basal_MIP.png';

imagesc(max(CiliaVolume(:,:,2,:),[],4));
colormap(green3.^0.5)
filename = '../Figures/Figure_1_Cilia_MIP.png';

imagesc(max(CiliaVolume(:,:,3,:),[],4));
colormap(blue3)
filename = '../Figures/Figure_1_DAPI_MIP.png';

imagesc(Output.FinalCombination_RGB)
filename = '../Figures/Figure_1_Output_RGB.png';
