afile = '..\data\annotation\blast_annotation\swiss_prot_go_annotation.txt';
[spe,id,~] = textread(afile,'%s%d%s');
length(unique(id))-5
specie_name = [];
for i=146011:length(spe)
spid = strsplit(spe{i},'_');
specie_name = [specie_name;spid(2)];
end