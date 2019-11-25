fpUTR = table;
%%
for i=1:24
    prefix = 'gene_TXCDUTR_ReadOutput_processedData_24filesDec_';
    fullname = [prefix, char(T.FileName(i)),'_shift15.txt'];
    file = readtable(fullname);
    if i==1
       fpUTR.AccNum = file.AccNum; 
    end
    setLabel = [char(T.Types(i)),'_',num2str(T.Timepoint(i))...
            ,'_',num2str(T.Sets(i))];
    fpUTR.(setLabel) = file.fpUTR_reads;
end
for i=25:30
    prefix = 'gene_TXCDUTR_ReadOutput_processedData_';
    fullname = [prefix, char(T.FileName(i)),'_shift15.txt'];
    file = readtable(fullname);
    if i==1
       fpUTR.AccNum = file.AccNum; 
    end
    setLabel = [char(T.Types(i)),'_',num2str(T.Timepoint(i))...
            ,'_',num2str(T.Sets(i))];
    fpUTR.(setLabel) = file.fpUTR_reads;
end

%%
fpUTR_w_RNA = fpUTR;
writetable(fpUTR_w_RNA,'fpUTR_with_RNA.csv');
%fpUTR = removevars(fpUTR, {'RNA_0_8','RNA_1_8','RNA_2_8','RNA_0_9','RNA_1_9','RNA_2_9'});

%% filter at 3

idx = min(fpUTR_w_RNA{:,2:end},[],2) >= 3 ; 
fpUTR_filtered = fpUTR_w_RNA(idx,:);

%% output
%writetable(fpUTR,'fpUTR_wo_RNA.csv');
writetable(fpUTR_filtered,'fpUTR_filtered3.csv');
