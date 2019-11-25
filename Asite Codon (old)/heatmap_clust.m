%% fake data
% try play with heatmap

load filteredyeastdata.mat
cgo = clustergram(yeastvalues(1:30,:),'Standardize','Row')
set(cgo,'RowLabels',genes(1:30),'ColumnLabels',times)
set(cgo,'Linkage','complete','Dendrogram',3)
set(cgo,'Colormap',redbluecmap)

cgo = clustergram(cod_comb{:,2:9},'standardize','Column')
set(cgo,'RowLabels',codon_names)
set(cgo,'ColumnLabels',cod_comb.Properties.VariableNames(2:end))


cgo = clustergram(cod_comb{:,2:end},...
                    'Cluster','All',...
                    'Standardize','Column',...
                    'RowLabels',codon_names,...
                    'Linkage','complete',...
                    'ColumnLabels',cod_comb.Properties.VariableNames(2:end))%,... 
                    %'RowLabelsColor',abc)
h = plot(cgo); set(h,'TickLabelInterpreter','none');
saveas(gcf,['codFreq_Asite_incl9/','Heatmap Codon Percent','.png']);