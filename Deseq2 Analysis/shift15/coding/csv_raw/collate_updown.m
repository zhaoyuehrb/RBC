S15ER2 = readtable('ER_S15_t2.csv');
up = S15ER2.log2FoldChange > 1;
down = S15ER2.log2FoldChange < -1;
writetable(S15ER2(up,1),'S15ER2_up.txt')
writetable(S15ER2(down,1),'S15ER2_down.txt')

