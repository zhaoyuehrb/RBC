% AccNum 2 GeneName
fgt = readtable('filteredGenesDetails_human_240118.txt');
acc2name = containers.Map;
for i = 1:height(fgt)
   acc2name(fgt.AccNum{i}) = fgt.GeneName{i}; 
end

dR = 'significant/';
files = dir(dR);
for i = 1:length(files)
    file = files(i).name;
    if length(file) < 5
       continue 
    end
    fT = readtable([dR,file],'Delimiter','\t','ReadVariableNames',0);
    names = table;
    for j = 1:height(fT)
       names.name{j} = acc2name(fT.Var1{j}); 
    end
    writetable(names,['geneName/',dR,file],'WriteVariableNames',0);
end