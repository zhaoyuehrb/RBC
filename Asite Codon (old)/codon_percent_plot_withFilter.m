if 0
    cod_comb = table;
    cod_comb.codon = codon_names;
    codon_pos = [16:18];
end

%filename = 'CGATGT-s_6_1 codon cut at 031.txt';
%curSet = 'S15_8';
t=readtable(filename);

counts_arr = zeros(64,1);
codon_map=containers.Map;
for i = 1:64
   codon_map(codon_names{i})=0;
end


tic 
seq_dat = t.Codon;
reads_dat = t.Reads;
acc_dat = t.AccNum;
skipped_count = 0;
unmapped_count = 0;
used_dat_count = 0;
maxLen = 0; % for debugging
for i = 1:length(reads_dat)
    read_raw = reads_dat(i);
    seq = seq_dat{i};
    codon = seq(codon_pos);
    accNum = acc_dat{i};
    % normalise by codon percentage && cdReads
    codon_idx = Codon2Idx(codon);
    accRow19210 = Acc2Row19210(accNum);
    codon_percent = codons_percent(accRow19210,codon_idx);
    % TODO: filter out wrong frame data
    lastLen = 18;
    
    cdORFLoc = strfind(ORFs(accNum),seq(1:lastLen));
    % multiple maps
    while length(cdORFLoc) > 1 & lastLen<29
       lastLen = lastLen+1;
       cdORFLoc = strfind(ORFs(accNum),seq(1:lastLen));
       if lastLen > maxLen
           maxLen = lastLen
       end
    end
    % unmapped
    if isempty(cdORFLoc)
        unmapped_count = unmapped_count+1;
        continue
    end
    % wrong frame as known
    if mod(cdORFLoc,3) ~= 1
        skipped_count = skipped_count+1;
        continue; 
    end
    % ignore wrong demand (shouldnt be here in any case
    % for debugging
    if codon_percent==0 
        %continue; 
        disp('aborted');
        return
    end
    cdRead = cdReads(Acc2Row8015(accNum));
    read_normalised = read_raw / (cdRead * codon_percent);
    if isKey(codon_map,codon)
       codon_map(codon) =  codon_map(codon) + read_normalised;
    else
        codon_map(codon) = read_normalised;
    end
    used_dat_count = used_dat_count +1;
    if mod(skipped_count,1000)==1
        disp(skipped_count)
    end
end
toc

for i = 1:64
   counts_arr(i)=codon_map(codon_names{i});
end
cod_comb.(curSet) = counts_arr;

figure('Renderer', 'painters', 'Position', [10 10 800 600])
[vals,idxs] = sort(counts_arr,'descend');
codon_ordered = codon_names(idxs);
x = categorical(codon_ordered);
x = reordercats(x,codon_ordered);
bar(x,vals/sum(vals));

title([curSet ' codon percentage'],'interpreter','none')
ylabel('Score')
xlabel('Codon')
%saveas(gcf,['codFreqIncl9/','Codon onRead Freq for ',curSet,'.png']);

disp(['Unmapped dat: ',int2str(unmapped_count)])
disp(['Used dat: ',int2str(used_dat_count)])
disp(['Total dat: ',int2str(length(reads_dat))])