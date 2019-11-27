
%bin_num = 'bin1';

for num = 1:5
bin_num = ['bin',num2str(num)];
label_flag = 0
AUGcodon = cell(24,1);
AUGbase = cell(24,1);
for i=16:18
    prefix = 'AUGbaseToReads_processedData_24filesDec_';
    fullname = [prefix, char(T.FileName(i)),'.txt'];
    AUGbase{i} = fullname;
    
    prefix = 'AUGcodonToReads_processedData_24filesDec_';
    fullname = [prefix, char(T.FileName(i)),'.txt'];
    AUGcodon{i} = fullname;
end

for i=25:27
    prefix = 'AUGcodonToReads_processedData_';
    fullname = [prefix, char(T.FileName(i)),'.txt'];
    AUGcodon{i-6} = fullname;
    
    prefix = 'AUGbaseToReads_processedData_';
    fullname = [prefix, char(T.FileName(i)),'.txt'];
    AUGbase{i-6} = fullname;
end

for i = [16:18,25:27]
    setLabel = [char(T.Types(i)),'_',num2str(T.Timepoint(i))...
            ,'_',num2str(T.Sets(i))];
    baseIdx = i;
    if i > 20
       baseIdx = i-6; 
    end
    basefile = readtable(['5bins/',bin_num,'/',AUGbase{baseIdx}]);
    codonfile = readtable(['5bins/',bin_num,'/',AUGcodon{baseIdx}]);
    figure;hold on;grid on;
    plot(basefile.Var1,basefile.Var2);
    
    if label_flag
        dat = basefile{1:end-1,:};
        offset = dat(1,1);
        % conditions
        [pks,locs] = findpeaks(dat(:,2),...
            'MinPeakDistance',6,'MinPeakHeight',max(dat(:,2))/8,'Threshold',0.08);

        text(locs+offset+.6,pks+.05,num2str(locs+offset))
        plot(locs+offset,pks,'o','MarkerSize',12)
    end
    ylabel('Read density (RPM)');
    xlabel('Distance from first start codon (nt)');
    title(['AUG base plot ',setLabel,' ',bin_num],'Interpreter','none')
    saveas(gcf,['5bins/plots/base/AUG_',setLabel,'_',bin_num,'.png']);
        
    
    figure;hold on;grid on;
    plot(codonfile.Var1,codonfile.Var2);
    
    if label_flag
        dat = codonfile{1:end-1,:};
        offset = dat(1,1);
        [pks,locs] = findpeaks(dat(:,2),...
            'MinPeakDistance',2,'MinPeakHeight',0.4);

        text(locs+offset+.6,pks+.05,num2str(locs+offset))
        plot(locs+offset,pks,'o','MarkerSize',12)
    end
    
    ylabel('Read density (RPM)');
    xlabel('Distance from first start codon (3nt)');
    title(['AUG codon plot ',setLabel,' ',bin_num],'Interpreter','none')
    
    saveas(gcf,['5bins/plots/codon/AUG_',setLabel,'_',bin_num,'.png']);
    
    
    
    
    close all
end  
end
