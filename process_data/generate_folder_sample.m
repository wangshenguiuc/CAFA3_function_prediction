specie_l = {'Human','Yeast','Drosophila','Mouse','Elegans'};
for sp =specie_l;
    specie = char(sp);
    [Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
    nnode = length(Gene_name);
    nfold = 3;
    rng(2)
nrepeat = 3;
    ntest = floor(nnode/3);
for i=1:nrepeat
    rp = randperm(nnode);
    st = 1;
    ed = ntest;
    test_ind = rp(st:ed);
    
    test_Gene_id = values(Gene_name_rev,num2cell(test_ind));
    
    fout = fopen(['../data/loo_gene_id/',specie,'_loo_test_gene_id_fold',num2str(i),'.txt'],'w');
    for g = test_Gene_id
        fprintf(fout,'%s\n',char(g));
    end
    fclose(fout);
end
end