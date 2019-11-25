if 1
codon_counts = containers.Map;
codon_bias   = containers.Map;

for i = 1:height(filterGeneTable)
    accNum = filterGeneTable.AccNum{i};
    codon_counts(accNum) = codoncount(ORFs(accNum));
    codon_bias(accNum)   = codonbias(ORFs(accNum));
end
end

%% access with codon_names{1}
codon_names = fieldnames(codon_counts(filterGeneTable.AccNum{1}));


%% Gene codon composition

codons_percent = zeros(height(filterGeneTable),64);
% stop codons at 49,51,57
for j = 1:64
    codon = codon_names{j};
    for i = 1:height(filterGeneTable)
        accNum          = filterGeneTable.AccNum{i};
        codon_len       = floor(ORF_len(accNum)/3);
        codon_freq      = codon_counts(accNum).(codon);
        if isnan(codon_freq)
            codons_percent(i,j) = 0;
        else
            codons_percent(i,j) = codon_freq/codon_len;
        end
    end
end


