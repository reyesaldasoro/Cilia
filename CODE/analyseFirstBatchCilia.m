
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

for k=1:numFiles
    tic
    disp(k)
    shortName{k}            = dir0(k).name(26:34);
    currFile                = strcat(baseDir,filesep,dir0(k).name);
    CiliaVolume             = readCilia(currFile);
    Output                  = segmentCilia(CiliaVolume,cp);
    % Save individual results
    Ratios_C_N(k,1)         = Output.Ratio_C_N;  
    Ratios_B_C(k,1)         = Output.Ratio_B_C;  
    q1                      = [Output.FinalCilia_MIP_P.MajorAxisLength];
    q2                      = q1(q1>0);
    LengthsPerCase(k,1:numel(q2))=q2;
    t2(k)=toc;
    %figure(k)
    %imagesc(Output.FinalCombination_RGB)
    % % Display results
    % h0=figure;
    % finalOutput(:,:,k)      = ((Output.FinalCilia_MIP==0).*Output.FinalNuclei_MIP)+(20+Output.FinalCilia_MIP);
    % h1=subplot(121);
    % imagesc(2*max(CiliaVolume(:,:,1:3,:)/16/255,[],4))
    % 
    % for k2=1:numel(Output.FinalCilia_MIP_P)
    %     currLength = Output.FinalCilia_MIP_P(k2).MajorAxisLength/ 4.8438;
    %     text(10+Output.FinalCilia_MIP_P(k2).Centroid(1),10+Output.FinalCilia_MIP_P(k2).Centroid(2),num2str(currLength,3),'color','w',FontSize=7)
    % end
    % title(strcat(shortName{k},',  ratio =',num2str(RatioPerCase(k))),'interpreter','none')
    % h2=subplot(122);
    % imagesc(finalOutput(:,:,k))
    % for k2=1:numel(Output.FinalCilia_MIP_P)
    %     currLength = Output.FinalCilia_MIP_P(k2).MajorAxisLength/ 4.8438;
    %     text(10+Output.FinalCilia_MIP_P(k2).Centroid(1),10+Output.FinalCilia_MIP_P(k2).Centroid(2),num2str(currLength,3),'color','w',FontSize=7)
    % end
    % colormap (jet2)
    % t2(k)=toc;
    % h0.Position = [ 488   309   829   353];
    % h1.Position=[0.05 0.06 0.44 0.88];
    % h2.Position=[0.55 0.06 0.44 0.88];
    % filename = strcat('Results/Res_2025_07_31_',shortName{k},'.png');
    % print('-dpng','-r100',filename)

end

    LengtsPerCaseC=LengthsPerCase/CalibrationFactor;


%%

h0=figure;
h1=gca;
hold off
plot(1:numFiles,Ratio,'r-o')
grid on

title('Maximum Intensities per Channel')
h1.TickLabelInterpreter = 'none';
h1.XTick                = 1:numFiles;
h1.XTickLabel           = shortName;
h1.XTickLabelRotation   = 90;
h1.Position             = [  0.0674    0.27    0.9112    0.60];
h0.Position             = [  643.4000   97.8000  812.0000  333.6000];




%%

figure
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
imagesc(CiliaVolume(:,:,2,4));colorbar

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



% ctRange = -5:1:6;
% 
% numResults = numel(ctRange);
% numRows = round(sqrt(numResults));
% 
% %%
% figure
% tiledlayout(numRows,ceil(numResults/numRows),TileSpacing="none",Padding="tight")
% 
% averageCellDiameter = 100;
% for ind = 1:numResults
%     disp(ind)
%     labels = segmentCells2D(cp,max(DAPI,[],3), ...
%         ImageCellDiameter=averageCellDiameter, ...
%         CellThreshold=ctRange(ind), ...
%         FlowErrorThreshold=10);        
% 
%     loverlay = labeloverlay(DAPI(:,:,k)/16,labels);
% 
%     nexttile
%     imshow(labels)    
%     title("Cell Threshold: " + num2str(ctRange(ind)))
% end
% linkaxes(findobj(gcf,Type="axes"))