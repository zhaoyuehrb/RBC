
codeMap = revgeneticcode;
aminos = fieldnames(codeMap);


cmap = lines(20);   % coloring
c_final = [ ];
x = [ ];
y = [ ];


for i = 2:21
    codes = (codeMap.(aminos{i}));
    for j = 1:length(codes)
        %codes{j}
        x=[x;[codes{j},'(',aminos{i},')']];
        y=[y,co_map(codes{j})];
        c_final = [c_final;cmap(i-1,:)];
    end
end
xcat = categorical(cellstr(x));
xcat = reordercats(xcat,cellstr(x));
b=bar(xcat,y);



b.FaceColor = 'flat';
% b.CData(2,:) = [.5 0 .5];
% for k = 1:length(y)
%     b.CData(k,:)= [k,k,k]/100;
% end
b.CData = c_final;
ax = gca;

for i = 1:61
   ax.XTickLabel{i} = sprintf('\\color[rgb]{%f,%f,%f}%s', c_final(i,:), ax.XTickLabel{i}); 
end

