addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
addpath '../data/'
addpath 'GeneMania'
specie_l = {'Elegans','Drosophila'};
for sp = specie_l
    specie = char(sp);
    [GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();

        [Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
        Gene_GO_train_annotation = read_annotation(specie,Gene_name,GO_name,GO_net);

    ngene = length(Gene_name);
    
    term_eia = pfp_eia(GO_net', logical(Gene_GO_train_annotation));
    %
    network = Gene_net;
    label_mat = Gene_GO_train_annotation;
    
    numNetworks = length(network);
    nnode =size(network{1},1);
    nlabel = size(label_mat,2);
    assert(nnode==size(label_mat,1));
    rng(2);
    % 7. get area under the AUC for prediction, as well as the scores
    pp = randperm(nnode); % this is the permutation index for doing the cross-validation;
    nFolds = 3; % number of cross-validation fold
    
    [areas, pr,weights,loo_sc,test_ind] = predictMClassesCV(label_mat, network, nFolds,pp);
    
    %
    % g1 = GO_namespace(:,1);
    % g2 = GO_namespace(:,2);
    % nlabel = size(label_mat,2);
    % go_type = unique(g2);
    % npos=zeros(1,nlabel);
    % st=[3,11,31,101,3];
    % ed=[10,30,100,300,300];
    % for i=1:nlabel
    %     npos(i)=nnz(label_mat(:,i));
    % end
    % for i=go_type'
    %     func_t=g1(g2==i);
    %     for j=1:5
    %         cat=func_t(npos(func_t)>=st(j) & npos(func_t)<=ed(j));
    %         fprintf('%d %d %d %f\n',i,j,length(cat),nanmean(nanmean(areas(i,cat))));
    %     end
    %
    % end
    
    eval_res = cell(1,nFolds);
        
    for p=1:1
        eval_res{p}.term_eval_res = term_evaluation(loo_sc(test_ind{p},:),label_mat(test_ind{p},:),label_mat,GO_namespace,0);
%         eval_res{p}.protein_eval_res = protein_evaluation(loo_sc(test_ind{p},:),label_mat(test_ind{p},:),term_eia,label_mat,GO_namespace);
    end
    save(['../inter_result/GeneMANIA/eval_res_select_code',specie]);
%    avg_mic_auroc = zeros(size(eval_res{1}.term_eval_res.mic_auroc));
%     avg_mac_auroc = zeros(size(eval_res{1}.term_eval_res.mac_auroc));
%     for p=1:length(mic_auroc)
%         avg_mic_auroc = avg_mic_auroc + eval_res{p}.term_eval_res.mic_auroc/length(eval_res{p}.term_eval_res.mic_auroc);
%         avg_mac_auroc = avg_mac_auroc + eval_res{p}.term_eval_res.mac_auroc/length(eval_res{p}.term_eval_res.mac_auroc);
%     end
%     dlmwrite(['..\result\auc\GM.',specie,'.auc'],[avg_mic_auroc;avg_mac_auroc]);
%     
end


