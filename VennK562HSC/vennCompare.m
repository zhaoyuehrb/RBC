pathK562 = 'Mapping/K562/deseq/shift15/csv_raw/fil2/';
pathHSC = 'deseq_yue/withRefSeq/shift15/coding/csv_raw/';
pathOut = 'VennK562HSC/';

fileK562 = 'res_S15_rpf.csv';
fileHSC = 'ER_S15_t2.csv';

tK562 = readtable([pathK562,fileK562]);
tHSC = readtable([pathHSC,fileHSC]);

%% direct compare

upsK562 = tK562.Var1(tK562.log2FoldChange>1);
upsHSC = tHSC.Var1(tHSC.log2FoldChange>1);

upsCommon = intersect(upsK562,upsHSC)

writetable(table(upsCommon),[pathOut,'K562 ',fileK562(5:end-4),...
    ' vs HSC ',fileHSC(1:end-4),' up all'],'WriteVariableNames',0);

downsK562 = tK562.Var1(tK562.log2FoldChange<-1);
downsHSC = tHSC.Var1(tHSC.log2FoldChange<-1);

downsCommon = intersect(downsK562,downsHSC)

writetable(table(downsCommon),[pathOut,'K562 ',fileK562(5:end-4),...
    ' vs HSC ',fileHSC(1:end-4), ' down all'],'WriteVariableNames',0);
%% compare within significant ones

upsK562 = tK562.Var1(tK562.log2FoldChange>1 & tK562.padj<0.05);
upsHSC = tHSC.Var1(tHSC.log2FoldChange>1 & tHSC.padj<0.05);

upsCommon = intersect(upsK562,upsHSC)

writetable(table(upsCommon),[pathOut,'K562 ',fileK562(5:end-4),...
    ' vs HSC ',fileHSC(1:end-4),' up sig'],'WriteVariableNames',0);

downsK562 = tK562.Var1(tK562.log2FoldChange<-1 & tK562.padj<0.05);
downsHSC = tHSC.Var1(tHSC.log2FoldChange<-1 & tHSC.padj<0.05);

downsCommonSig = intersect(downsK562,downsHSC)

writetable(table(downsCommonSig),[pathOut,'K562 ',fileK562(5:end-4),...
    ' vs HSC ',fileHSC(1:end-4), ' down sig'],'WriteVariableNames',0);






