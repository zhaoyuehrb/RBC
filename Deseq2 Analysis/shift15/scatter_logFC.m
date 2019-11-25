fpNames = dir('shift15/csv_Raw/fil10/*.csv');
for i = 1:length(fpNames)
    fpName = fpNames(i).name;
    fpFile = readtable(['shift15 fpUTR/csv_raw/',fpName]);
    cdFile = readtable(['shift15/csv_raw/fil10/',fpName]);
    [commonGene,i_fp,i_cd] = intersect(fpFile.Var1,cdFile.Var1);
    lFC_fp = fpFile.log2FoldChange(i_fp);
    lFC_cd = cdFile.log2FoldChange(i_cd);
    uL = max([lFC_fp;lFC_cd]);
    lL = min([lFC_fp;lFC_cd]);
    
    % filters for points to label
    idx_label = lFC_fp > prctile(lFC_fp,99.9) | ...
                lFC_fp < prctile(lFC_fp,0.1) | ...
                lFC_cd > prctile(lFC_cd,99.9) | ...
                lFC_cd < prctile(lFC_cd,0.1);
            
    idxs = find(idx_label==1);

    figure;hold on;grid on;
    scatter(lFC_fp,lFC_cd);
    plot([lL,uL],[lL,uL],'LineStyle',':');
    xlabel('fpUTR log2FoldChange')
    ylabel('exon log2FoldChange')
    for j = idxs'
         text(lFC_fp(j),lFC_cd(j)-rand()*1,...
             [' ',acc2name(commonGene{j})],'FontSize',12)
    end

    ylim([lL,uL]);
    xlim([lL,uL]);
    title([fpName(1:end-4),' Scatter Plot'],'Interpreter','none')
    saveas(gcf,['compareBoth/',fpName(1:end-4),' Scatter Plot.jpg']);

end

    
%     x=lFC_fp(idxs);
%     y=lFC_cd(idxs) - rand(length(idxs),1);
%     labels = {};
%     for j = 1:length(idxs)
%     labels{j} = acc2name(commonGene{idxs(j)});
%     end
%     labelpoints(x, y, labels, 'N', 1, 1)




