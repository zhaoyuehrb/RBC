if 0
for i = 2:31
    RPKMs{strcmp(RPKMs{:,i},'-'),i} = {'0'};
end


for i=2:31
    RPKMs.(i) = str2double(RPKMs{:,i});
end
end

%% TEs
if 0    % get logFC from deseq csv results, muz order in RPKMs sequence
        %                               same as filterGeneTable sequence
    csv_path = '../deseq_yue/Proper Version/csv_raw/';
    deseq_table = readtable('../deseq_yue/Proper Version/csv_raw/TE_t2.csv');
    gene_list = deseq_table.Var1;
    filter_idx = ismember(RPKMs.AccNum,gene_list);
    gene_list_ordered = RPKMs.AccNum(filter_idx);
    oidxs = zeros(length(gene_list),1);
    for i = 1:length(gene_list)
        oidxs(i) = find(strcmp(gene_list,gene_list_ordered{i}));
    end
    TE = deseq_table.log2FoldChange(oidxs);
    cur_set = 'TE t2';
end

if 1
% filter_idx = RPKMs.RNA_0_9 >= 10 & RPKMs.totalRPF_2_8 >= 10;
% TE = RPKMs.S24_2_8(filter_idx) ./ RPKMs.totalRPF_2_8(filter_idx);
% cur_set = 'RNA t0 set9';
 
 filter_idx = RPKMs.RNA_2_9 >= 10;
 TE = RPKMs.RNA_2_9(filter_idx);
 cur_set = 'RNA t2 set9';


% filter_idx = RPKMs.totalRPF_0_8 >= 10 & RPKMs.totalRPF_2_8 >= 10 &...
%             RPKMs.RNA_0_8 >= 10 & RPKMs.RNA_2_8 >= 10;
% 
% TE_0 = RPKMs.totalRPF_0_8(filter_idx) ./ RPKMs.RNA_0_8(filter_idx);
% TE_1 = RPKMs.totalRPF_2_8(filter_idx) ./ RPKMs.RNA_2_8(filter_idx);
% TE = log(TE_1 ./ TE_0);
%      
% cur_set = 'TE t2 vs t0';
end

co = zeros(64,1);
co_map = containers.Map;
for i=1:64
    if 0
        figure;scatter(codons_percent(filter_idx,i),TE_total);
        title([cur_set,' codon ', codon_names{i}]);
        xlabel('codon percentage');
        ylabel('corr');
        saveas(gcf,['totalRPF_codon_scatter/',cur_set,' codon ', codon_names{i},'.png']);
        close all
    end
    A = codons_percent(filter_idx,i);
    r = corrcoef(A,TE);
    co(i) = r(2,1);
    co_map(codon_names{i}) = r(2,1);
end

%plot(co)

if 1    % bar plot ordered by correlation
    codon_wo_stop_idx = [1:48,50,52:56,58:64];
    [vals,idxs] = sort(co(codon_wo_stop_idx),'descend');
    original_idx = codon_wo_stop_idx(idxs);
    x = categorical(codon_names(original_idx));
    x = reordercats(x,codon_names(original_idx));
    bar(x,vals);
    title([cur_set,' corr against codon sorted']);
    ylabel('corr coef');
    xlabel('codons');
    saveas(gcf,['Codon Corr Sorted/',cur_set,'.png']);
end

if 0  % bar plot ordered by amino acid
    codeMap = revgeneticcode;
    aminos = fieldnames(codeMap);
    cmap = lines(20);   % coloring
    c_final = [ ];
    x = [ ];
    y = [ ];
    for i = 2:21
        codes = (codeMap.(aminos{i}));
        for j = 1:length(codes)
            %codes{j}
            x=[x;[codes{j},'(',aminos{i},')']];
            y=[y,co_map(codes{j})];
            c_final = [c_final;cmap(i-1,:)];
        end
    end
    xcat = categorical(cellstr(x));
    xcat = reordercats(xcat,cellstr(x));
    
    figure('Renderer', 'painters', 'Position', [10 10 1200 800])
    b=bar(xcat,y);
    b.FaceColor = 'flat';
    % b.CData(2,:) = [.5 0 .5];
    % for k = 1:length(y)
    %     b.CData(k,:)= [k,k,k]/100;
    % end
    b.CData = c_final;
    ax = gca;
    for i = 1:61
       ax.XTickLabel{i} = sprintf('\\color[rgb]{%f,%f,%f}%s', c_final(i,:), ax.XTickLabel{i}); 
    end
    title(cur_set);
    ylabel('corr coef');
    xlabel('codons/amino acid');
    %saveas(gcf,['TEamino/',cur_set,'.png']);
end


if 1  % bar plot ordered by amino acid, sorted based corr within each amino
    codeMap = revgeneticcode;
    aminos = fieldnames(codeMap);
    cmap = lines(20);   % coloring
    c_final = [ ];
    x = [ ];
    y = [ ];
    for i = 2:21
        codes = (codeMap.(aminos{i}));
        x_tmp = [];
        y_tmp = [];
        for j = 1:length(codes)
            %codes{j}
            % add sort, idea: store in tmp, sort, then add to x and y
            x_tmp = [x_tmp;[codes{j},'(',aminos{i},')']];
            y_tmp = [y_tmp,co_map(codes{j})];
            c_final = [c_final;cmap(i-1,:)];
        end
        if length(codes) > 1
            [y_tmp,sort_idx] = sort(y_tmp,'descend');
             x_tmp = x_tmp(sort_idx,:);
        end
        x=[x;x_tmp];
        y=[y,y_tmp];
    end
    xcat = categorical(cellstr(x));
    xcat = reordercats(xcat,cellstr(x));
    
    figure('Renderer', 'painters', 'Position', [10 10 1200 800])
    b=bar(xcat,y);
    b.FaceColor = 'flat';
    % b.CData(2,:) = [.5 0 .5];
    % for k = 1:length(y)
    %     b.CData(k,:)= [k,k,k]/100;
    % end
    b.CData = c_final;
    ax = gca;
    for i = 1:61
       ax.XTickLabel{i} = sprintf('\\color[rgb]{%f,%f,%f}%s', c_final(i,:), ax.XTickLabel{i}); 
    end
    title([cur_set,' corr against codon by amino acid']);
    ylabel('corr coef');
    xlabel('codons/amino acid');
    saveas(gcf,['Codon Corr by Amino Acid/',cur_set,'.png']);
end
%close all
