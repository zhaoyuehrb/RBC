FASTAData = fastaread('splicedORFs_240118.txt');
ORFs = containers.Map;
for i=1:length(FASTAData)
    header = FASTAData(i).Header;
    newStr = split(header,'_chr');
    acc_num = char(newStr(1));
    if isKey(ORFs,acc_num)
       disp(acc_num) 
    end
    ORFs(acc_num) = FASTAData(i).Sequence;
end


filterGeneTable = readtable('filteredGenesDetails_human_240118.txt');
ORF_len = containers.Map;
ORF_len_list = zeros(19210,1);

for i = 1:height(filterGeneTable)
    accNum = filterGeneTable.AccNum{i};

    curr_len = filterGeneTable.ORFLength(i);
    if curr_len ~= length(ORFs(accNum))
        disp([accNum,' ',curr_len,'  ',length(ORFs(accNum))]);
    end
    ORF_len(accNum) = curr_len;
    ORF_len_list(i) = curr_len;
end