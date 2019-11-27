
for num = 5:5

label_flag = 0
AUGcodon = cell(24,1);
AUGbase = cell(24,1);
for i=18:18
    prefix = 'AUGbaseToReads_processedData_24filesDec_';
    fullname = [prefix, char(T.FileName(i)),'.txt'];
    AUGbase{i} = fullname;
    
    prefix = 'AUGcodonToReads_processedData_24filesDec_';
    fullname = [prefix, char(T.FileName(i)),'.txt'];
    AUGcodon{i} = fullname;
end

for i=27:27
    prefix = 'AUGcodonToReads_processedData_';
    fullname = [prefix, char(T.FileName(i)),'.txt'];
    AUGcodon{i-6} = fullname;
    
    prefix = 'AUGbaseToReads_processedData_';
    fullname = [prefix, char(T.FileName(i)),'.txt'];
    AUGbase{i-6} = fullname;
end

for i = [18:18,27:27]
    setLabel = [char(T.Types(i)),'_',num2str(T.Timepoint(i))...
            ,'_',num2str(T.Sets(i))];
    baseIdx = i;
    if i > 20
       baseIdx = i-6; 
    end
    basefile = readtable(['t2normed/thr50/',AUGbase{baseIdx}]);
    codonfile = readtable(['t2normed/thr50/',AUGcodon{baseIdx}]);
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
    ylabel('Read density');
    xlabel('Distance from first start codon (nt)');
    title(['AUG base plot ',setLabel],'Interpreter','none')
    saveas(gcf,['t2normed/thr50/AUG_',setLabel,'_base_normed.png']);
        
    
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
    
    ylabel('Read density');
    xlabel('Distance from first start codon (3nt)');
    title(['AUG codon plot ',setLabel],'Interpreter','none')
    
    saveas(gcf,['t2normed/thr50/AUG_',setLabel,'_codon_normed.png']);
    
    
    
    
    close all
end  
end