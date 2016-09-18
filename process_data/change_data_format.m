mat = dlmread('..\data\embedding_vector\clusDCA\Yeast_Gene_dim500 rsp_0.5.US');
specie = 'Yeast'
% [Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);

[Trans_ID,~] = textread('..\analysis\build_Meng_input\input\_Y.ppi.node.txt','%s%s');

Trans_ID_num = cell2mat(values(Gene_name,Trans_ID));


dlmwrite('..\data\embedding_vector\clusDCA\Yeast_Gene_500_init',mat(Trans_ID_num,:)','delimiter','\t')
