addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
specie = 'Human';
use_clusDCA_embedding = false;
read_exist_embedding = true;
use_Meng_embedding = ~use_clusDCA_embedding;
local_test = false;
result_file = '../result/function_prediction/Human_model.txt';
auc_result_file = '../result/function_prediction/Human_model_detail_AUC.txt';
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
 
[read_file_list,s1,s2,s3,s4,s5,s6] = textread(result_file,'%s%f%f%f%f%f%f','delimiter','\t');
select_file = find(s1>0.89);
read_file_list = read_file_list(select_file);
F = dir('..\data\embedding_vector\LINE\9_type_vec\*.*');
for ii = 1:length(F)
    if isempty(find(strcmp(read_file_list,F(ii).name), 1))
        continue;
    end
    fname = F(ii).name;
    embedding_file_name =fname;
    Gene_embedding = read_embedding(Gene_name, embedding_file_name);    
    model_name = [fname]
    [eval_res,mic_auroc,mac_auroc]=clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation, term_eia,GO_namespace,specie,GO_net);
%     fprintf('bilinear %s mic = %f, mac=%f\n',F(ii).name,mean(mic_auroc),mean(mac_auroc));
    write_result_to_file(model_name,mic_auroc(:),mac_auroc(:),auc_result_file);
end

% 
% GO_embedding_file_name = '../data/embedding_vector/clusDCA/GO.embedding';
% GO_embedding = read_embedding_GO( GO_name, GO_embedding_file_name);
% [read_file_list,~,~,~,~,~,~] = textread(result_file,'%s%f%f%f%f%f%f','delimiter','\t');
% F = dir('..\data\embedding_vector\LINE\9_type_vec\*.*');
% for ii = 79:length(F)
%     if ~isempty(find(strcmp(read_file_list,F(ii).name), 1))
%         continue;
%     end
%     if ~isempty(find(strcmp(F(ii).name,'_enum'), 1))
%         continue;
%     end
%     fname = F(ii).name;
%     fname_l = strsplit(fname,'_');
%     mode = str2num(fname_l{length(fname_l)});
%     if isempty(mode)
%         continue
%     end
%     fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_GO_M.ant_Y.ant_D.ant_E.ant.node_500_10_0.800000_3250000_ppi_rwr_0_seq_rwr_0_enum_9_lr_0.025000_train_mode_2';
%     embedding_file_name =fname;
%     Gene_embedding = read_embedding(Gene_name, embedding_file_name);    
%     model_name = [fname];
%     [eval_res,mic_auroc,mac_auroc]=clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation, term_eia,GO_namespace,specie,GO_net);
%     fprintf('bilinear %s mic = %f, mac=%f\n',F(ii).name,mean(mic_auroc),mean(mac_auroc));
%     write_result_to_file(model_name,mic_auroc,mac_auroc,result_file);
% end
% 


% model_name = 'embedding seq kNN';
% [mic_auroc,mac_auroc]=blast(embed_net,Gene_GO_annotation(:,:), GO_namespace);
% write_result_to_file(model_name,mic_auroc,mac_auroc,result_file);
%
% model_name = 'blast';
% [mic_auroc,mac_auroc]=blast(Seq_net,Gene_GO_annotation, GO_namespace);
% write_result_to_file(model_name,mic_auroc,mac_auroc,result_file);

%
% if use_Meng_embedding
%     if read_exist_embedding
%         model_name = 'CAFA2';
%         embedding_file_name = 'Human_Mouse_Yeast_Drosophila_Elegans_seq_Mouse.annot_Yeast.annot_Drosophila.annot_Elegans.annot.node_500_10_0.800000_300'
%         %  embedding_file_name = '../data/embedding_vector/Meng_vec1_6_species/vec500.emb';
%         Gene_embedding = read_embedding(Gene_name, embedding_file_name,1,' ');
%     else
%         run_Meng_embedding(Gene_name_rev,Gene_net,specie);
%     end
% else
%     if read_exist_embedding
%         model_name = 'clusDCA 500';
%         embedding_file_name = '../data/embedding_vector/clusDCA/Gene_dim500 rsp_0.5.US'
%         Gene_embedding = dlmread(embedding_file_name);
%     else
%         Gene_embedding = learn_DCA_vector(Gene_net,0.5,500,Gene_name,specie);
%     end
% end
% if read_exist_embedding
%     %     GO_embedding_file_name = '../data/embedding_vector/clusDCA/hp.embedding'
%     GO_embedding_file_name = '../data/embedding_vector/clusDCA/GO.embedding';
%     %     GO_embedding = dlmread(GO_embedding_file_name);
%     GO_embedding = read_embedding_GO( GO_name, GO_embedding_file_name);
%     %     GO_embedding = learn_DCA_vector(GO_net,0.8,1500,GO_name,'GO');
% end
% if local_test
%     model_name = 'embedding seq bilinear multi species';
%     [mic_auroc,mac_auroc]=clusDCA(Gene_embedding, GO_embedding, Gene_GO_annotation, GO_namespace,specie);
%     write_result_to_file(model_name,mic_auroc,mac_auroc,result_file);
%     read_result_from_file
% else
%     [Gene_GO_train_annotation,Gene_GO_test_annotation]  = read_CAFA2_annotation(specie,Gene_name,GO_name,GO_net);
%     test_clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation,Gene_GO_test_annotation, GO_namespace,specie);
% end