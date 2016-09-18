addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
specie = 'Drosophila'; %HMYDE
% specie_prefix = 'H.ant_M.ant_Y.ant_E.ant.node.txt';

[GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();

[Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
Gene_GO_train_annotation = read_annotation(specie,Gene_name,GO_name,GO_net);

ngene = length(Gene_name);

term_eia = pfp_eia(GO_net, logical(Gene_GO_train_annotation));


GO_embedding_file_name = '../data/embedding_vector/clusDCA/GO_dim2500 rsp_0.8.US';
GO_embedding = dlmread(GO_embedding_file_name);

    
result_file = ['../result/function_prediction/',specie,'_model_selected.txt'];
fid = fopen(result_file, 'a+');
fclose(fid);
[read_file_list,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16] = textread(result_file,'%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter','\t');
F = dir(['..\data\embedding_vector\selected_code\',specie,'\*.*']);
for ii = 1:length(F)
    if length(F(ii).name) <3
        continue
    end
    if ~isempty(find(strcmp(read_file_list,F(ii).name), 1))
        continue;
    end
    fname = F(ii).name
    embedding_file_name =['\\?\E:\swang141\project\SequencingNetwork\Sheng\data\embedding_vector\selected_code\',specie,'\',fname];
    Gene_embedding = read_embedding(Gene_name, embedding_file_name);
    model_name = [fname];
    [eval_res,mic_auroc,mac_auroc,std_mac_auroc]=clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation, term_eia,GO_namespace,specie,GO_net,10);
    %     fprintf('bilinear %s mic = %f, mac=%f\n',F(ii).name,mean(mic_auroc),mean(mac_auroc));
    write_result_to_file(model_name,mic_auroc,mac_auroc,result_file);

% write_auc_result_to_file([specie,'_',method],mic_auroc,mac_auroc,std_mac_auroc);
end
