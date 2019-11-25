files = dir('../csv_filtered/*.csv');
x = NaN(3,10000);   % rows are  tputr(deplete,enrich),...
                    %           fp
                    %           orf
                    % cols are lengths
                    % for boxplot
j=0;
for file = files'
    rowIdx=mod(j,3)+1;
    t=readtable(['../csv_filtered/',file.name]);
    fp_table = t(:,1);
    for i = 1:height(t)
        str_orf = ORFs(char(fp_table.Var1(i)));
        codon_orf = codoncount(str_orf);
        % taken from https://proteinstructures.com/Structure/Structure/amino-acids.html
        hydroPAminos = ['A';'I';'L';'M';'F';'V';'P';'G'];
        hydros = 0;
        for ami = hydroPAminos'
            for don = 1:length(codeMap.(ami))
                hydros = hydros + codon_orf.(codeMap.(ami){don});
            end
        end


        x(rowIdx,i) = hydros * 3 / length(str_orf);
    end
    if rowIdx==3
        medians=nanmedian(x');
        x_max=nanmax(x');
        figure;hold on;
        H=boxplot(x','Labels',{'ORF down','ORF NA','ORF up'});
        x=NaN(3,10000);
        get(H,'tag');
        set(H(6,:),'color','k','linewidth',0.5)
        %color=['r','g','b','k','y','m','c','r','g'];
        color='g';
        h=findobj(gca,'Tag','Box');
        for i=1:length(h)
            patch(get(h(i),'XData'),get(h(i),'YData'),color,'FaceAlpha',0.5);
        end

        for i=1:3
            text(i-0.1,x_max(i)+0.02,num2str(medians(i)));
        end
        title([file.name(1:end-7),'_hydrophobicity_boxplot'],'Interpreter', 'none');
        saveas(gcf,['hydrophobicity/',file.name(1:end-7),'_hydrophobicity_boxplot.png']);
    end
    
    j=j+1;
end
