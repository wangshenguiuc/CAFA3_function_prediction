addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
specie_l = {'Human','Yeast','Elegans','Drosophila','Mouse'};
specie_nick_name_l = {'Human','Yeasx','CAEEL','Drome','Rat'};
taxid_l = {'9606','4932','6239','7227','10090'};
sp_ind = 1;
specie = specie_l{sp_ind};
specie_nick_name = specie_nick_name_l{sp_ind};
taxid = taxid_l{sp_ind};
[GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
[Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
Gene_GO_train_annotation = read_annotation(specie,Gene_name,GO_name,GO_net);
term_eia = pfp_eia(GO_net, logical(Gene_GO_train_annotation));



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

method ='Blast_clusDCA_combined_prediction';
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
final_score_our =zscore(blast_final_score);
[mic_auroc,mac_auroc,homo_mac_auroc_detail,std_mac_auroc] = evaluation(final_score_our(test_ind,:),Gene_GO_train_annotation(test_ind,:),Gene_GO_train_annotation,GO_namespace,1);

[~,~,~,~,clusDCA_final_score_base]=clusDCA(Gene_embedding_base, GO_embedding, Gene_GO_train_annotation, term_eia,GO_namespace,specie,GO_net,-1);
final_score_base = zscore(clusDCA_final_score_base);
[mic_auroc,mac_auroc,ppi_mac_auroc_detail,std_mac_auroc] = evaluation(final_score_base(test_ind,:),Gene_GO_train_annotation(test_ind,:),Gene_GO_train_annotation,GO_namespace,1);

final_score_integ = zscore(final_score_base)+zscore(final_score_our);
[mic_auroc,mac_auroc,integ_mac_auroc_detail,std_mac_auroc] = evaluation(final_score_integ(test_ind,:),Gene_GO_train_annotation(test_ind,:),Gene_GO_train_annotation,GO_namespace,1);

for thres = 0.8:0.05:0.95
for i=1:2
    number = zeros(4,4);

           ntotal_all = 0;
    acc_ppi_all = 0;
    acc_homo_all = 0;
    acc_integ_all = 0;
    both_ppi_homo_all = 0;
    for j=1:4         
        ntotal = length(ppi_mac_auroc_detail{i,j});
        acc_ppi =length(find(ppi_mac_auroc_detail{i,j}>thres));
        acc_homo =length(find(homo_mac_auroc_detail{i,j}>thres));
        acc_integ = length(find(integ_mac_auroc_detail{i,j}>thres));
        both_ppi_homo = length(intersect(find(ppi_mac_auroc_detail{i,j}>thres),find(homo_mac_auroc_detail{i,j}>thres)));
        number(j,:) = [both_ppi_homo/ntotal,(acc_homo)/ntotal,(acc_ppi)/ntotal,(acc_integ)/ntotal];
        ntotal_all = ntotal_all + ntotal;
        acc_ppi_all = acc_ppi_all + acc_ppi;
        acc_homo_all = acc_homo_all + acc_homo;
        acc_integ_all = acc_integ_all  + acc_integ;
        both_ppi_homo_all = both_ppi_homo_all + both_ppi_homo;
 fprintf('thres:%f %d %d ppi:%f homo:%f both:%f integ:%f\n',thres,i,j,(acc_ppi)/ntotal,(acc_homo)/ntotal,both_ppi_homo/ntotal,(acc_integ)/ntotal);

    end
       fprintf('thres:%f %d ppi:%f homo:%f both:%f integ:%f\n',thres,i,acc_ppi_all,acc_homo_all,both_ppi_homo_all,acc_integ_all);
 
    legend = {'Intersection of network and homology','Homology','Network','Integrated'};
    grouplabel = {'3-10','11-30','31-100','101-300'};
    if i==1
        title = [specie,' MF'];
        name = [specie,'_MF_comple_',num2str(thres*100)];
    else
        title = [specie,' BP'];
        name = [specie,'_BP_comple_',num2str(thres*100)];
    end
    xlabel = 'Number of annoated genes';
    ylabel = ['Percentage of functions with AUROC > ',num2str(thres)];
    plot_bar_complementary(title,legend,number,xlabel,ylabel,grouplabel,name,(i-1)*2+2);
end
end