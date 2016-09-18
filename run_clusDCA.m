addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'

specie_l = {'Human'};
for sp = specie_l
    specie = char(sp);
    use_clusDCA_embedding = false;
    read_exist_embedding = true;
    use_Meng_embedding = ~use_clusDCA_embedding;
    local_test = false;
    result_file = '../result/function_prediction/Human_model.txt';
    [GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
    if strcmp(specie,'all')
        [Seq_net,Gene_name,Gene_name_rev] = read_sequence_network();
        Gene_GO_annotation = read_annotation(specie,Gene_name,GO_name,GO_net);
    else
        [Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
        %     Seq_net= read_sequence_network(Gene_name);
        Gene_GO_train_annotation = read_annotation(specie,Gene_name,GO_name,GO_net);
    end
    ngene = length(Gene_name);
    
    term_eia = pfp_eia(GO_net, logical(Gene_GO_train_annotation));
    
     GO_embedding_file_name = '../data/embedding_vector/clusDCA/GO_dim2500 rsp_0.8.US';   
     GO_embedding = dlmread(GO_embedding_file_name);%learn_DCA_vector(GO_net,0.8,2500,GO_name,specie);
    % embedding_file_name = '../data/embedding_vector/clusDCA/Gene_dim2500 rsp_0.5.US';
    Gene_embedding = learn_DCA_vector(Gene_net,0.5,2500,Gene_name,specie);
    
%     [eval_res,mic_auroc,mac_auroc]=clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation, term_eia,GO_namespace,specie);
%     avg_mic_auroc = zeros(size(mic_auroc{1}));
%     avg_mac_auroc = zeros(size(mac_auroc{1}));
%     for p=1:length(mic_auroc)
%         avg_mic_auroc = avg_mic_auroc + mic_auroc{p}/length(mic_auroc{p});
%         avg_mac_auroc = avg_mac_auroc + mac_auroc{p}/length(mac_auroc{p});
%     end
%     dlmwrite(['..\result\auc\clusDCA.',specie,'.auc'],[avg_mic_auroc;avg_mac_auroc]);
end