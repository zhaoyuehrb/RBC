%% complete table for cdReads>=10
cdReadsTable = readtable('Fil10_BioRep_ReadCounts.txt');


%% first AccNum 2 row number map
AccNums = cdReadsTable.AccNum;
Acc2Row19210 = containers.Map;
Acc2Row8015 = containers.Map;
for i = 1: length(AccNums)
    Acc2Row19210(AccNums{i}) = find(strcmp(filterGeneTable.AccNum,AccNums{i}));
    Acc2Row8015(AccNums{i}) = i;
end

%% normalise cdReads and codon percentage
% actually done alr,
% just access using codon_percent(Acc2Row(accNum))
% codons are under codon_names{123}
% for many accesses of codon
Codon2Idx = containers.Map;
for i = 1:length(codon_names)
   Codon2Idx(codon_names{i}) = i;
end
%% sel 1 cdReads for demo
cdReads = cdReadsTable.Total0cdRead_set9;
