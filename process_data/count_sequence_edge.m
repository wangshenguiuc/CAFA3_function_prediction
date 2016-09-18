specie_l = {'Human','Yeast','Drosophila','Mouse','Elegans'};
[g1,g2,~,~] = textread('../analysis/process_sequence_network/sequence_network.txt','%s%s%s%s');
g1=  upper(g1);
g2 = upper(g2);
for sp1=1:5
    specie = char(specie_l(sp1));
    [Gene_name,~,~] = read_string_network(specie);
    gene_set{sp1} = keys(Gene_name);
end

for sp1=1:5
    [ia1,~] = ismember(g1,gene_set{sp1});
    specie1 =  char(specie_l(sp1));
    for sp2 = 1:5
        specie2 =  char(specie_l(sp2));
        [ia2,~] = ismember(g2,gene_set{sp2});
        ia = ia1 & ia2;
        nseq = length(find(ia));
        fprintf('%s %s %d\n',specie1,specie2,nseq);
    end
end