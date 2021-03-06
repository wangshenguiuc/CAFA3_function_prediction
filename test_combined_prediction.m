addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
specie_l = {'Human','Yeast','Elegans','Drosophila','Mouse'};
specie_nick_name_l = {'Human','Yeasx','CAEEL','Drome','Rat'};
taxid_l = {'9606','4932','6239','7227','10090'};
for sp_ind =1:5
    specie = specie_l{sp_ind};
    specie_nick_name = specie_nick_name_l{sp_ind};
    taxid = taxid_l{sp_ind};
    [GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
    [Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
    Gene_GO_train_annotation = read_annotation(specie,Gene_name,GO_name,GO_net);
    
    GO_embedding_file_name = '../data/embedding_vector/clusDCA/GO_dim2500 rsp_0.8.US';
    GO_embedding = dlmread(GO_embedding_file_name);%learn_DCA_vector(GO_net,0.8,2500,GO_name,specie);
    %     embedding_file_name = ['../data/embedding_vector/clusDCA/',specie,'_Gene_dim2500 rsp_0.5.US'];
    %     fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_GOloo_H_M.asc_Y.asc_D.asc_E.asc.node.txt_500_10_0.500000_50_ppi_rwr_1_seq_rwr_1_enum_9_lr_0.025000_mode_1_meta__H7M9G_H7E9G_H7Y9G_H7D9G_H9G8G_H9G';
    %
    %     our_embedding_file_name = ['\\?\E:\swang141\project\SequencingNetwork\Sheng\data\embedding_vector\selected_code\',specie,'\',fname];
    %     Gene_embedding_our = read_embedding(Gene_name, our_embedding_file_name);
    %
    base_embedding_file_name =[specie,'_Gene_dim2500 rsp_0.5.US'];
    Gene_embedding_base = dlmread(base_embedding_file_name);
    
    method ='Blast_clusDCA_combined_prediction_strong';
    [nnode,nlabel] = size(Gene_GO_train_annotation);
    nfold = 3;
    ntest = nnode/nfold;
    rng(2);
    rp = randperm(nnode);
    st = 1;
    ed = ntest;
    test_ind = rp(st:ed);
    train_ind = rp(ed+1:nnode);
    
    Gene_name_mul_spe = cell(1,5);
    for sp=1:5
        spt = char(specie_l(sp));
        [Gene_name_mul_spe{sp},~] = read_string_network(spt);
    end
    
    [seq_net,Gene_name_new,Gene_name_new_rev] = read_sequence_network_subnet(Gene_name_mul_spe,test_ind,sp_ind);
    clear Gene_name_mul_spe;
    All_Gene_GO_annotation = read_annotation('all',Gene_name_new,GO_name,GO_net);
    
    blast_score =  blast(seq_net, All_Gene_GO_annotation);
    blast_final_score = zeros(nnode,nlabel);
    blast_final_score(test_ind,:) = blast_score;
    % [~,~,~,~,clusDCA_final_score_our]=clusDCA(Gene_embedding_our, GO_embedding, Gene_GO_train_annotation, term_eia,GO_namespace,specie,GO_net,10);
    final_score_blast =zscore(blast_final_score);
    [mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc] = evaluation(final_score_blast(test_ind,:),Gene_GO_train_annotation(test_ind,:),Gene_GO_train_annotation,GO_namespace,1);
    
    method = 'blast';
    write_auc_result_to_file([specie,'_',method],mic_auroc,mac_auroc,std_mac_auroc,specie);
    %PSB submision step = 50
    [~,~,~,~,clusDCA_final_score_base]=clusDCA(Gene_embedding_base, GO_embedding, Gene_GO_train_annotation, term_eia,GO_namespace,specie,GO_net,10);
    final_score_clusDCA = zscore(clusDCA_final_score_base);
    [mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc] = evaluation(final_score_clusDCA(test_ind,:),Gene_GO_train_annotation(test_ind,:),Gene_GO_train_annotation,GO_namespace,1);
    
    method = 'clusDCA_opt';
    write_auc_result_to_file([specie,'_',method],mic_auroc,mac_auroc,std_mac_auroc,specie);
    
    %PSB submision step = -1
    
    %     blast_score_norm = (blast_final_score-min(blast_final_score(:)))/(max(blast_final_score(:))-min(blast_final_score(:)));
    %     embed_score_norm = (clusDCA_final_score_base-min(clusDCA_final_score_base(:)))/(max(clusDCA_final_score_base(:))-min(clusDCA_final_score_base(:)));
    %
    final_score_add =zscore(final_score_blast)+zscore(final_score_clusDCA);
    [mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc] = evaluation(final_score_add(test_ind,:),Gene_GO_train_annotation(test_ind,:),Gene_GO_train_annotation,GO_namespace,1);
    method = 'blast_clusDCA_zscore';
    write_auc_result_to_file([specie,'_',method],mic_auroc,mac_auroc,std_mac_auroc,specie);
    
    final_score_add =zscore(final_score_blast)+zscore(final_score_clusDCA);
    [mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc] = evaluation(final_score_add(test_ind,:),Gene_GO_train_annotation(test_ind,:),Gene_GO_train_annotation,GO_namespace,1);
    method = 'blast_clusDCA_zscore_trans';
    write_auc_result_to_file([specie,'_',method],mic_auroc,mac_auroc,std_mac_auroc,specie);
end