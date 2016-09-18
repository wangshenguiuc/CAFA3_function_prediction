specie_l={'Human','Mouse','Yeast','Drosophila','Elegans'};
for species = specie_l
    specie = char(species);
    [Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
    
    nnet = length(Gene_net);
    network_list = {'coexpression','cooccurrance','database','experimental','fusion','neighborhood'};
    
    for i=1:nnet
        adj_network = Gene_net{i};
        Gene_net{i} = diag(sum(adj_network,2).^-0.5)*adj_network*diag(sum(adj_network,2).^-0.5);
        fout = fopen(['../data/network/STRING/',specie,'.',network_list{i},'.dn'],'w');
        [a,b,c] = find(Gene_net{i});
        prefix = ['STRING_',network_list{i}];
        for j=1:length(c)
            fprintf(fout,'%s\t%s\t%f\t%s\n',char(values(Gene_name_rev,{a(j)})),char(values(Gene_name_rev,{b(j)})),...
                c(j),prefix);
        end
        fclose(fout);
    end
end

