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
        x(rowIdx,i) = length(fpUTRs(char(fp_table.Var1(i))));
        x(rowIdx+3,i) = length(tpUTRs(char(tp_table.Var1(i))));
        x(rowIdx+6,i) = length(ORFs(char(tp_table.Var1(i))));
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
            text(i-0.3,x_max(i)+1000,num2str(medians(i)));
        end
        title([file.name(1:end-7),' nt length_boxplot'],'Interpreter', 'none');
        saveas(gcf,['lengths/',file.name(1:end-7),' UTR length boxplot.png']);
    end
    
    j=j+1;
end
