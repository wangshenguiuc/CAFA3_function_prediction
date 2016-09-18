addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
specie = 'Human';
% for go_type = {'bp','cc','mf','hp'}
GO_type = 'bp';
[GO_name, GO_net, GO_name_rev] = read_CAFA2_GO_network(GO_type);
% % [CC_GO_name, CC_GO_net] = read_CAFA2_GO_network('cc');
% [MF_GO_name, MF_GO_net] = read_CAFA2_GO_network('mf');
% % [HP_GO_name, HP_GO_net] = read_CAFA2_GO_network('hp');
% % 
% % model_name = 'CAFA2';
embedding_file_name = 'CAFA2_target_t1_seq.node.txt_300_10_0.800000_500';
[Gene_embedding,Gene_name, Gene_name_rev] = read_embedding_new_Gene(embedding_file_name);
% 
% 
GO_embedding_file_name = ['../data/embedding_vector/clusDCA/',GO_type,'.embedding'];
GO_embedding = read_embedding_GO(GO_name,GO_embedding_file_name);
% GO_embedding = learn_DCA_vector(GO_net,0.8,500,GO_name,go_type);
% end
% 
[Gene_GO_train_annotation,Target_gene,Gene_embedding,Gene_name,Gene_name_rev]  = read_CAFA2_annotation(Gene_name,GO_name,GO_net,Gene_embedding,Gene_name_rev);
% mac_auroc = cv_clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation);
% mac_auroc = predict_clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation,Target_gene)
% 
final_score = test_clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation, Target_gene,'kNN');
% [mic_auroc,mac_auroc]=clusDCA(Gene_embedding, GO_embedding, Gene_GO_annotation, GO_namespace,specie);
% final_score = sum(Gene_GO_train_annotation)/max(sum(Gene_GO_train_annotation));
CAFA2_save_data(final_score,GO_type);
CAFA2_evaluate;
load('E:\swang141\project\SequencingNetwork\Sheng\analysis\CAFA2\github\CAFA2-master\evaluation\hpo_all_type1_mode1\hpo_all_type1_mode1\our_most_popular.mat')
seq_fmax
mean(term_auc.auc(~isnan(term_auc.auc)))