specie_l = {'Human','Yeast','Drosophila','Mouse','Elegans'};
specie_nick_name_l = {'HUMAN','YEAST','DROME','RAT','CAEEL'};
taxnomy_id = {'9606','4932','7227','10090','6239'};

fout =fopen('../data/sequence_data/other_uniprot_species_to_5species.txt','w');
for sp=1:5
    specie = char(specie_l(sp));
    specie_nick_name = char(specie_nick_name_l(sp));
    tid = char(taxnomy_id(sp));
    [Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
    gid = keys(Gene_name);
    gid = upper(gid);
    [ensp,alias,~] = textread(['../data/sequence_data/',char(tid),'.protein.aliases.v10.txt'],'%s%s%s','headerlines',1);
    alias = upper(alias);
    
    [LIA,LOC] = ismember(gid,alias);
    
    selected_ensp = ensp(LOC(LOC~=0));
    selected_ensg = gid(LIA);
    
    ensp_to_ensg = containers.Map(selected_ensp,selected_ensg);
    
    [ensp,uniprot,evalue,~,~,~] = textread(['../data/sequence_data/',char(tid),'_to_uniprot.evalue.txt'],'%s%s%f%f%f%f','headerlines',1);
    evalue(evalue < 1e-100) = 1e-100;
    filter = evalue<1e-8 & isKey(ensp_to_ensg,ensp);
    ensp = ensp(filter);
    evalue = evalue(filter);
    uniprot = uniprot(filter);
    filter = zeros(length(uniprot),1);
    for x = specie_nick_name_l
        cx = char(x);
        filter = filter + (1-cellfun(@isempty,strfind(uniprot,cx)));
    end
    filter = filter==0 & isKey(ensp_to_ensg,ensp);
    ensg =values(ensp_to_ensg, ensp(filter));
    uniprot = uniprot(filter);
    logevalue = -1*log10(evalue(filter))/100;
    for i=1:length(logevalue)
        fprintf(fout,'%s\t%s\t%f\n',char(ensg(i)),char(uniprot(i)),logevalue(i));
    end
end
fclose(fout);

