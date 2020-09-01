T = sortrows(T,[3,4]);

%% add file name column to T

rpkmFilName = cell(30,1);
for i=1:30
    prefix = 'shift15_geneTXCD_RPKMoutput_';
    fullname = [prefix, char(T.FileName(i)),'_shift15.txt'];
    rpkmFilName{i} = fullname;
end

T.rpkmFilName = rpkmFilName;
%% populate 2 tables

cdReadsTable = table();
cdRPKMTable = table();
file = readtable(T.rpkmFilName{1});
cdReadsTable.AccNum = file.AccNum;
cdRPKMTable.AccNum = file.AccNum;
for i = 1:30
    file = readtable(T.rpkmFilName{i});
    setLabel = [char(T.Types(i)),'_',num2str(T.Timepoint(i))...
            ,'_',num2str(T.Sets(i))];
    cdReadsTable.(setLabel) = file.cdReads;
    cdRPKMTable.(setLabel) = file.cdRPKM;
end

%% format rpkm table (actually no need.. sua


%% cdReads filter

idxs = sum((cdReadsTable{:,2:end}>=5)') == 30;
disp(['num of genes pass filter: ' num2str(sum(idxs))])

cdReadsTableFiltered = cdReadsTable(idxs,:);
%writetable(cdReadsTableFiltered,'cdReadsFiltered8.csv')

%writetable(cdReadsTableFiltered,'cdReadsNotFiltered.csv');

cdRPKMTableFiltered = cdRPKMTable(idxs,:);
%writetable(cdRPKMTableFiltered,'cdRPKMTableFiltered8.csv')

rpkm = table();
rpkm.AccNum = cdRPKMTable.AccNum;
cols = cdRPKMTable.Properties.VariableNames;
for i = 2:2:30
    colname = cols{i}(1:end-2)
    rpkm.(colname) = (cdRPKMTable{:,i} + cdRPKMTable{:,i+1} ) / 2;
end
writetable(rpkm,'cdRPKMCombFil5.csv')