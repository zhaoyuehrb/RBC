rpkm_over_cd = zeros(19210,30);
%%
for i=1:24
    filename = char(T.FileName(i));
    rpkms = readtable(['../coding/shift15_geneTXCD_RPKMoutput_24filesDec_',filename,'_shift15.txt']);
    % reorder rpkms accnums
    idxs = zeros(19210,1);
    for j=1:19210
        idxs(j)=find(strcmp(rpkms.AccNum,fpUTR.AccNum{j}));
    end
    rpkms = rpkms(idxs,:);    
    rpkms.cdRPKM(strcmp(rpkms.cdRPKM,'-')) = {'0'};
    rpkms.cdRPKM = str2double(rpkms.cdRPKM);
    rpkm_over_cd(:,i) = rpkms.cdRPKM ./ rpkms.cdReads;
end
for i=25:30
    filename = char(T.FileName(i));
    rpkms = readtable(['../coding/shift15_geneTXCD_RPKMoutput_processedData_',filename,'_shift15.txt']);
    % reorder rpkms accnums
    idxs = zeros(19210,1);
    for j=1:19210
        idxs(j)=find(strcmp(rpkms.AccNum,fpUTR.AccNum{j}));
    end
    rpkms = rpkms(idxs,:);    
    rpkms.cdRPKM(strcmp(rpkms.cdRPKM,'-')) = {'0'};
    rpkms.cdRPKM = str2double(rpkms.cdRPKM);
    rpkm_over_cd(:,i) = rpkms.cdRPKM ./ rpkms.cdReads;
end


%%
%fpUTR_w_RNA = fpUTR;
%fpUTR = removevars(fpUTR, {'RNA_0_8','RNA_1_8','RNA_2_8','RNA_0_9','RNA_1_9','RNA_2_9'});

fpUTR_rpkm = fpUTR_w_RNA;
for i = 1:30
    fpUTR_rpkm{:,i+1} = fpUTR_rpkm{:,i+1} .* rpkm_over_cd(:,i);
end

writetable(fpUTR_rpkm,'fpUTR_all_rpkm.csv');