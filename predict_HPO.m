addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
specie = 'hpo';
GO_type = 'hpo';
result_file = '../result/function_prediction/HPO_model.txt';
[GO_name, GO_name_rev,GO_net, GO_namespace] = read_HPO_network();
[Gene_name,Gene_name_rev,Gene_net] = read_string_network('human');

Gene_GO_train_annotation = read_annotation(specie,Gene_name,GO_name,GO_net,'hpoa.annotation.t0');
embedding_file_name = 'Human_Mouse_Yeast_Drosophila_Elegans_seq_Mouse.annot_Yeast.annot_Drosophila.annot_Elegans.annot.node_500_20_0.800000_500';
Gene_embedding = read_embedding(Gene_name, embedding_file_name);

GO_embedding_file_name = '../data/embedding_vector/clusDCA/hp.embedding';
GO_embedding = read_embedding_GO( GO_name, GO_embedding_file_name);

Gene_GO_test_annotation = read_annotation(specie,Gene_name,GO_name,GO_net,'hpoa.annotation.t1');
Gene_GO_test_annotation = Gene_GO_test_annotation - Gene_GO_train_annotation;
Gene_GO_test_annotation = double(Gene_GO_test_annotation>0);

[human_target,CAFA2_all_target_id,CAFA2_hpo_target_id]= textread('../data/CAFA2_target_gene/hpo_cafa2_target_ensg.txt','%s%d%d');
human_target_id = cell2mat(values(Gene_name,human_target));

final_score = test_clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation,Gene_GO_test_annotation, GO_namespace,specie,human_target_id);
for i=1:size(final_score,2)
    final_score(:,i) = final_score(:,i)-mean(final_score(:,i));
end
pred_obj = CAFA2_save_data(final_score,GO_type,CAFA2_hpo_target_id);
CAFA2_evaluate;
load('E:\swang141\project\SequencingNetwork\Sheng\analysis\CAFA2\github\CAFA2-master\evaluation\hpo_all_type1_mode1\hpo_all_type1_mode1\our_most_popular.mat')
seq_fmax
mean(term_auc.auc(~isnan(term_auc.auc)))
our_pred= load(['..\analysis\CAFA2\github\CAFA2-master\baselines\hpo\our_most_popular'],'pred');

our_eval_cafa=load('E:\swang141\project\SequencingNetwork\Sheng\analysis\CAFA2\github\CAFA2-master\evaluation\hpo_all_type1_mode1\hpo_all_type1_mode1\our_most_popular.mat');
our_oa=load(['..\analysis\CAFA2\github\CAFA2-master\benchmark\groundtruth\bpoa.mat']);
% mac_auroc=mac_auc_evaluation(pred_obj.score(CAFA2_hpo_target_id,:), our_oa.annotation(CAFA2_all_target_id,:))
% mac_auroc=mac_auc_evaluation(pred.score(CAFA2_hpo_target_id,:), oa.annotation(CAFA2_all_target_id,:))
[Lia,Locb] = ismember(our_pred.pred.object,our_oa.oa.object);

our_eva = protein_evaluation(our_pred.pred.score,our_oa.oa.annotation(Locb,:));
