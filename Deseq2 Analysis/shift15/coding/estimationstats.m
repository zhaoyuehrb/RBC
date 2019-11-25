dat=readtable('csv_raw/ER_S24_t0.csv');

if 0    % load gene list from csv (in this case Erythroid)
    gene_list = readtable('Gene Groups/Erythroid_06032019List.txt');
    filterGeneDetail = readtable('filteredGenesDetails_human_240118.txt');
    acc_list = {};

    for i = 1:height(gene_list)
       gene_name = upper(gene_list.Gene{i});
       idx = find(strcmp(filterGeneDetail.GeneName,gene_name));
       if ~isempty(idx)
            acc_list{end+1}=filterGeneDetail.AccNum{idx};
       end
    end
end

idx = ismember(dat.Var1,acc_list); %selected groups
identifiers = cell(8015,1);
identifiers(idx) = {'Erythroids'};
identifiers(~idx) = {'non-Erythroids'};
T = table(identifiers,dat.log2FoldChange,'VariableNames',{'Groups','logFC'});
writetable(T,'csv_dabest/ER_S24_t0_erythroid.csv')
%dabest('csv_dabest/ER_S24_t0_erythroid.csv')
