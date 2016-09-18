function  [hold_set,qseqid,valid_pos] = find_hold_out_gene_name(test_gene_id,oa_gene_id,specie_id)

oa_gene_id = upper(oa_gene_id);
test_gene_id = upper(test_gene_id);
[ensp,alias,~] = textread(['../data/sequence_data/',specie_id,'.protein.aliases.v10.txt'],'%s%s%s','headerlines',1);
alias = upper(alias);

[valid_pos,LOC] = ismember(test_gene_id,alias);

qseqid = ensp(LOC(LOC~=0));

[~,LOC] = ismember(oa_gene_id,alias);
oa_ensp = ensp(LOC(LOC~=0));
gene_map = containers.Map(oa_ensp, alias(LOC(LOC~=0)));

hold_out_ensp = intersect(oa_ensp,qseqid);

hold_out_gene = values(gene_map,hold_out_ensp);

[hold_set,~] = ismember(oa_gene_id,hold_out_gene);


valid_gene = values(gene_map,hold_out_ensp);

[hold_set,~] = ismember(oa_gene_id,hold_out_gene);