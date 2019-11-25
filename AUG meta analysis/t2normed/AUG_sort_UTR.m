t = readtable('t2normed/peakUTR_geneList_GGCTAC-s_7_1.txt')
[res,idx] = sort(t.Var2);

AccNums = t.Var1(idx(19210:-1:19190));

for i = 1:length(AccNums)
    geneNames{i} = acc2name(AccNums{i});
end

topGenes = table();
topGenes.AccNum = AccNums;
topGenes.Symbol = geneNames';
writetable(topGenes,'topGenesAtUTR.csv')

