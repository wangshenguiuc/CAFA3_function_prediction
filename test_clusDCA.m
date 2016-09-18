addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
method = 'clusDCA_old_GOA';
specie_l = {'Mouse'};
for i=1:1
specie = specie_l{i};
[GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
[Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
Gene_GO_train_annotation = read_annotation(specie,Gene_name,GO_name,GO_net);
ngene = length(Gene_name);

term_eia = pfp_eia(GO_net, logical(Gene_GO_train_annotation));

GO_embedding_file_name = '../data/embedding_vector/clusDCA/GO_dim2500 rsp_0.8.US';
GO_embedding = dlmread(GO_embedding_file_name);%learn_DCA_vector(GO_net,0.8,2500,GO_name,specie);
embedding_file_name = ['../data/embedding_vector/clusDCA/',specie,'_Gene_dim2500 rsp_0.5.US'];
Gene_embedding = dlmread(embedding_file_name);

[eval_res,avg_mic_auroc,avg_mac_auroc]=clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation, term_eia,GO_namespace,specie,GO_net,-1);

write_auc_result_to_file([specie,'_',method],avg_mic_auroc,avg_mac_auroc);
end