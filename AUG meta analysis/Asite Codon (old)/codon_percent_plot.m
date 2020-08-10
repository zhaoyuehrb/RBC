%filename = 'CGATGT-s_7_1 codon cut at 031.txt';
%t=readtable(filename);
counts_arr = zeros(64,1);
codon_map=containers.Map;
for i = 1:64
   codon_map(codon_names{i})=0;
end

tic 
codons_dat = t.Codon;
reads_dat = t.Reads;
acc_dat = t.AccNum;
for i = 1:length(reads_dat)
    read_raw = reads_dat(i);
    codon = codons_dat{i}(16:18);
    accNum = acc_dat{i};
    % normalise by codon percentage && cdReads
    codon_idx = Codon2Idx(codon);
    %accRow19210 = Acc2Row19210(accNum);
    %codon_percent = codons_percent(accRow19210,codon_idx);
    % ignore wrong demand
%     if codon_percent==0 
%         continue; 
%     end
    % TODO: filter out wrong frame data
    
    
    
    
    
    %cdRead = cdReads(Acc2Row8015(accNum));
    %read_normalised = read_raw / (cdRead * codon_percent);
    if isKey(codon_map,codon)
       codon_map(codon) =  codon_map(codon) + read_raw;
    else
        codon_map(codon) = read_raw;
    end
end
toc

figure('Renderer', 'painters', 'Position', [10 10 800 600])
codons = keys(codon_map);
for i = 1:64
   counts_arr(i)=codon_map(codons{i});
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