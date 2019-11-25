files = dir('../csv_filtered/*.csv');
x = NaN(9,10000);   % rows are  tputr(deplete,enrich),...
                    %           fp
                    %           orf
                    % cols are lengths
                    % for boxplot
j=0;
for file = files'
    rowIdx=mod(j,3)+1;
    t=readtable(['../csv_filtered/',file.name]);
    fp_table = t(:,1);
    tp_table = fp_table;
    for i = 1:height(t)
        str_tmp = fpUTRs(char(fp_table.Var1(i)));
        
        x(rowIdx,i) = (length(strfind(str_tmp,'G')) + length(strfind(str_tmp,'C')))/length(str_tmp);
        str_tmp = tpUTRs(char(fp_table.Var1(i)));
        x(rowIdx+3,i) = (length(strfind(str_tmp,'G')) + length(strfind(str_tmp,'C')))/length(str_tmp);
        str_tmp = ORFs(char(fp_table.Var1(i)));
        x(rowIdx+6,i) = (length(strfind(str_tmp,'G')) + length(strfind(str_tmp,'C')))/length(str_tmp);
    end
    if rowIdx==3
        medians=nanmedian(x');
        x_max=nanmax(x');
        figure;hold on;
        H=boxplot(x','Labels',{'fp down','fp NA','fp up',...
                'tp down','tp NA','tp up',...
                'ORF down','ORF NA','ORF up'});
        x=NaN(9,10000);
        get(H,'tag');
        set(H(6,:),'color','k','linewidth',0.5)
        %color=['r','g','b','k','y','m','c','r','g'];
        color='g';
        h=findobj(gca,'Tag','Box');
        for i=1:length(h)
            patch(get(h(i),'XData'),get(h(i),'YData'),color,'FaceAlpha',0.5);
        end

        for i=1:9
            text(i-0.3,x_max(i)+0.02,num2str(medians(i)));
        end
        title([file.name(1:end-7),'_GC_content_boxplot'],'Interpreter', 'none');
        saveas(gcf,['GCcontent/',file.name(1:end-7),'_GC_content_boxplot.png']);
    end
    
    j=j+1;
end
