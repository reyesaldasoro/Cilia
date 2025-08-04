for k1=1:numFiles
    h(k1) = subplot(6,6,k1);
    imagesc((basalbodies(:,:,k1)));
    axis off
    h(k1).Position(4)=0.13;
    h(k1).Position(3)=0.13;
end

%%

for k1=1:numFiles
    imagesc((basalbodies(:,:,k1)));colorbar
    filename = strcat('Results/BasalBody_',num2str(k1),'_MIP.png');
    title(shortName{k1},'interpreter','none');
        print('-dpng','-r100',filename)

end