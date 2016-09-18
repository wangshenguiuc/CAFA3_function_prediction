specie_l = {'Human','Yeast','Drosophila','Mouse','Elegans'};
specie_nick_name_l = {'HUMAN','YEAST','DROME','RAT','CAEEL'};
taxnomy_id = {'9606','4932','7227','10090','6239'};

for sp=1:5
    specie = char(specie_l(sp));
    specie_nick_name = char(specie_nick_name_l(sp));
    tid = char(taxnomy_id(sp));
    
    [Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
    [GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
    Gene_GO_train_annotation = read_annotation(specie,Gene_name,GO_name,GO_net);
    
    fout =fopen(['../data/annotation/',char(specie),'.GO.complete'],'w');
    [gid,goid,s] = find(Gene_GO_train_annotation);
    for i=1:length(s)
        fprintf(fout,'%s\t%s\t1\tGO_term\n',char(values(GO_name_rev,num2cell(goid(i)))),char(values(Gene_name_rev,num2cell(gid(i)))));
    end
    fclose(fout);
end

afile = '..\data\annotation\blast_annotation\swiss_prot_go_annotation.txt';
oa = read_blast_annotation(afile,GO_name,GO_net,0,specie_nick_name);

filter = zeros(1,length(oa.object));
 for x = specie_nick_name_l
        cx = char(x);
        filter = filter + (1-cellfun(@isempty,strfind(oa.object,cx)));
    end
filter = filter==0;

oa.object(filter) = [];
oa.annotation(filter,:) = [];

Gene_GO_train_annotation = oa.annotation;

fout =fopen('../data/annotation/blast_annotation/other_uniprot_species_to_5species.GO','w');
[gid,goid,s] = find(Gene_GO_train_annotation);
for i=1:length(s)
         fprintf(fout,'%s\t%s\t1\tGO_term\n',char(oa.ontology.term(goid(i))),char(oa.object(gid(i))));
 end
fclose(fout);

