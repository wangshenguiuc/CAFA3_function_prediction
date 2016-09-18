addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
specie = 'Elegans';
method = 'GeneMANIA.select_code';
load(['../inter_result/GeneMANIA/eval_res_select_code',specie]);
eval_res = cell(1,nFolds);
avg_mic_auroc = [];
avg_mac_auroc = [];
for p=1:1
    [mic_auroc,mac_auroc,~,std_mac_auroc] = evaluation(loo_sc(test_ind{p},:),label_mat(test_ind{p},:),label_mat,GO_namespace,0);
    %         eval_res{p}.protein_eval_res = protein_evaluation(loo_sc(test_ind{p},:),label_mat(test_ind{p},:),term_eia,label_mat,GO_namespace);
avg_mic_auroc = [avg_mic_auroc, mic_auroc];
avg_mac_auroc = [avg_mac_auroc, mac_auroc];
end

write_auc_result_to_file([specie,'_',method],avg_mic_auroc,avg_mac_auroc,std_mac_auroc,specie);
