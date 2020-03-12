filename = 'CGATGT-s_7_1 codon cut at 1821.txt';
t=readtable(filename);

tic 
codon_counts=containers.Map;
codons_dat = t.Codon;
reads_dat = t.Reads;
for i = 1:length(reads_dat)
    
    if isKey(codon_counts,codons_dat{i})
       codon_counts(codons_dat{i}) =  codon_counts(codons_dat{i}) + reads_dat(i);
    else
        codon_counts(codons_dat{i}) = reads_dat(i);
    end
end
toc

figure('Renderer', 'painters', 'Position', [10 10 800 600])
codons = keys(codon_counts);
for i = 1:64
   counts_arr(i)=codon_counts(codons{i});
end

[vals,idxs] = sort(counts_arr,'descend');
codon_ordered = codons(idxs);
x = categorical(codon_ordered);
x = reordercats(x,codon_ordered);
bar(x,vals/sum(vals));

title(filename)
ylabel('Freq')
xlabel('Codon')
%saveas(gcf,['Included Plots/','Codon onRead Freq for Total RPF set9 T0','.png']);