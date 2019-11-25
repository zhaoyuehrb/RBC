cod_comb = table;
cod_comb.codon = codon_names;
codon_pos = 13:15;
for ii = [16:18,25:30]
    filename = ['codFreqIncl9/',T.FileName{ii},' codon cut at 031.txt'];
    curSet = [T.Types{ii},'_',num2str(T.Timepoint(ii))...
        ,'_',num2str(T.Sets(ii))];
    codon_percent_plot_withFilter;
end


if 0

cod2Ami = geneticcode;

% color same amino
ami2color = containers.Map;
cmap = lines(21);
%cmap = cmap(1:3:63,:);
c_idx = 1;
for i = 1:64
   ami = cod2Ami.(codon_names{i});
   if ~isKey(ami2color,ami)
       ami2color(ami) = cmap(c_idx,:);
       c_idx = c_idx+1;
   end
end

% add amino acid name
c_final = zeros(64,3);
for i = 1:64    
    xlabels{i} = [codon_names{i} ,'(', cod2Ami.(codon_names{i}),')'];
    c_final(i,:) = ami2color(cod2Ami.(codon_names{i}));
end

for i = 2:25
counts_arr = cod_comb{:,i};
curSet = cod_comb.Properties.VariableNames{i}
figure('Renderer', 'painters', 'Position', [10 10 800 600])

[vals,idxs] = sort(counts_arr,'descend');
codon_ordered = xlabels(idxs);
c_ordered = c_final(idxs,:);
x = categorical(codon_ordered);
x = reordercats(x,codon_ordered);
b = bar(x,vals/sum(vals));

%b.FaceColor = 'flat';
%b.CData = c_ordered;
ax = gca;
for ii = 1:64
   ax.XTickLabel{ii} = sprintf('\\color[rgb]{%f,%f,%f}%s', c_ordered(ii,:), ax.XTickLabel{ii}); 
end


title([curSet ' codon percentage'],'interpreter','none')
xlabel('Codon')
ylabel('Score')


saveas(gcf,['codFreq_Asite_incl9/','Codon Freq for ',curSet,'.png']);
end
end