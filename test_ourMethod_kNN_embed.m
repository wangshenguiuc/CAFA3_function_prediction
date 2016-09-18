addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
specie_l = {'Human','Yeast','Drosophila','Mouse','Elegans'};
test_sid = 3;
specie = specie_l{test_sid};
[GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
Gene_name = cell(1,5);
Gene_name_rev = cell(1,5);
Gene_net = cell(1,5);
Gene_GO_train_annotation = cell(1,5);
Gene_embedding = cell(1,5);

GO_embedding_file_name = '../data/embedding_vector/clusDCA/GO_dim2500 rsp_0.8.US';
GO_embedding = dlmread(GO_embedding_file_name);
% fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_GOloo_Y_H.asc_M.asc_D.asc_E.asc.node.txt_500_10_0.500_170_1_1_9_0.010_2_1__Y7E9G_Y7H9G_Y7M9G_Y7D9G_Y9G8G_Y9G';
%fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_GOloo_M_H.asc_Y.asc_D.asc_E.asc.node.txt_500_10_0.500_7_1_1_9_0.025_2_1__M7E9G_M7H9G_M7Y9G_M7D9G_M9G8G_M9G';
% fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_GOloo_H_M.asc_Y.asc_D.asc_E.asc.node.txt_500_10_0.500000_50_ppi_rwr_1_seq_rwr_1_enum_9_lr_0.025000_mode_1_meta__H7M9G_H7E9G_H7Y9G_H7D9G_H9G8G_H9G';
% fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_GOloo_E_H.asc_M.asc_Y.asc_D.asc.node.txt_500_10_0.500_10_1_1_9_0.010_2_1__E7M9G_E7H9G_E7Y9G_E7D9G_E9G8G_E9G';
fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_GOloo_D_H.asc_M.asc_Y.asc_E.asc.node.txt_500_10_0.500_20_1_1_9_0.010_2_1__D7M9G_D7H9G_D7Y9G_D7E9G_D9G8G_D9G'
embedding_file_name =['\\?\E:\swang141\project\SequencingNetwork\Sheng\data\embedding_vector\selected_code\',specie,'\',fname];

for sp=1:5
    spt = char(specie_l(sp));
    [Gene_name{sp},Gene_name_rev{sp}] = read_string_network(spt);
    Gene_GO_train_annotation{sp} = read_annotation(spt,Gene_name{sp},GO_name,GO_net);
    Gene_embedding{sp} = read_embedding(Gene_name{sp}, embedding_file_name);
end

base_embedding_file_name =  ['../data/embedding_vector/clusDCA/',specie,'_Gene_dim2500 rsp_0.5.US'];
Gene_embedding_base = dlmread(base_embedding_file_name);
other_gene_go_annotation = [];
other_embedding = [];
for sp=1:5
    if sp==test_sid
        continue
    end
    other_gene_go_annotation = [other_gene_go_annotation;Gene_GO_train_annotation{sp}];
    other_embedding = [other_embedding;Gene_embedding{sp}];
end

lx = Gene_embedding{test_sid};
nspecie = length(Gene_embedding);
nlabel = size(GO_embedding,1);
nnode = cell(1,nspecie);
for i=1:nspecie
    nnode{i}= size(Gene_embedding{i},1);
end

nfold = 3;

rng(2)

ntest = floor(nnode{test_sid}/nfold);
nrepeat = 1;
st_d = 1:10:501;
ed_d = 11:10:511;
npos=zeros(1,nlabel);
rp = randperm(nnode{test_sid});
st = 1;
ed = ntest;
test_ind = rp(st:ed);
train_ind = rp(ed+1:nnode{test_sid});

test_embedding = lx(test_ind,:);
same_spe_embedding = lx;
other_spe_embedding = other_embedding;
fprintf('start to get distnace matrix');
same_D = pdist2(test_embedding,same_spe_embedding, 'cosine');
other_D = pdist2(test_embedding,other_spe_embedding, 'cosine');

same_ant = Gene_GO_train_annotation{test_sid};
same_ant(test_ind,:) = 0;
other_score=zeros(nnode{test_sid},nlabel);
final_score=zeros(nnode{test_sid},nlabel);
same_score=zeros(nnode{test_sid},nlabel);

other_ant = other_gene_go_annotation;
all_D = [same_D,other_D];
all_ant = [same_ant;other_ant];
all_score=zeros(nnode{test_sid},nlabel);

% base_score=zeros(nnode{test_sid},nlabel);
% [base_D,Gene_name_new,Gene_name_new_rev] = read_sequence_network_subnet(Gene_name,test_ind,test_sid);
% base_D=1-base_D;
% All_Gene_GO_annotation = read_annotation('all',Gene_name_new,GO_name,GO_net);
% base_ant=All_Gene_GO_annotation;
st=[3,11,31,101,3];
ed=[10,30,100,300,300];

npos = zeros(1,nlabel);
nvote = 4000;
for i=1:nlabel
    npos(i)=nnz(Gene_GO_train_annotation{test_sid}(:,i));
end
[~,~,~,~,embed_score_base]=clusDCA(Gene_embedding_base, GO_embedding, same_ant, [],GO_namespace,specie,GO_net,10);

[~,~,~,~,embed_score]=clusDCA(Gene_embedding{test_sid}, GO_embedding, same_ant, [],GO_namespace,specie,GO_net,10);
nvote = 4000;
    for category=1:2
        func_t1=GO_namespace(GO_namespace(:,2)==category,1);
            for i = 1:ntest
                [v, o] = sort(same_D(i, train_ind));
                o = o(~isinf(v));
                k = min(nvote, length(o));
                votes = sum(bsxfun(@rdivide, same_ant(train_ind(o(1:k)),func_t1)', same_D(i, train_ind(o(1:k)))), 2);
                same_score(test_ind(i),func_t1) = votes;
            end
            for i = 1:ntest
                [v, o] = sort(other_D(i, train_ind));
                o = o(~isinf(v));
                k = min(nvote, length(o));
                votes = sum(bsxfun(@rdivide, other_ant(train_ind(o(1:k)),func_t1)', other_D(i, train_ind(o(1:k)))), 2);
                other_score(test_ind(i),func_t1) = votes;
            end
            %         for i = 1:ntest
            %             [v, o] = sort(all_D(i, train_ind));
            %             o = o(~isinf(v));
            %             k = min(nvote, length(o));
            %             votes = sum(bsxfun(@rdivide, all_ant(train_ind(o(1:k)),func_t)', all_D(i, train_ind(o(1:k)))), 2);
            %             all_score(test_ind(i),func_t) = votes;
            %         end
    
    end
    fprintf('nvote=%d\n',nvote);
    
    final_score = zscore(same_score)+zscore(other_score)+zscore(embed_score)+zscore(embed_score_base);
    fprintf('start same\n')
    [mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc]  = evaluation(final_score(test_ind,:),Gene_GO_train_annotation{test_sid}(test_ind,:),Gene_GO_train_annotation{test_sid},GO_namespace,1);
    method = 'our_method';
    write_auc_result_to_file([specie,'_',method],mic_auroc,mac_auroc,std_mac_auroc,specie);
    
%     final_score = zscore(embed_score_base)+zscore(other_score);
%     fprintf('start both\n')
%     [mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc]  = evaluation(final_score(test_ind,:),Gene_GO_train_annotation{test_sid}(test_ind,:),Gene_GO_train_annotation{test_sid},GO_namespace,1);
% same_score_norm = (same_score-min(same_score(:)))/(max(same_score(:))-min(same_score(:)));
% other_score_norm = (other_score-min(other_score(:)))/(max(other_score(:))-min(other_score(:)));
% embed_score_norm = (embed_score-min(embed_score(:)))/(max(embed_score(:))-min(embed_score(:)));
% final_score = zscore(same_score_norm+other_score_norm+embed_score_norm)+zscore(same_score)+zscore(other_score)+zscore(embed_score);
%     [mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc]  = evaluation(final_score(test_ind,:),Gene_GO_train_annotation{test_sid}(test_ind,:),Gene_GO_train_annotation{test_sid},GO_namespace,1);
%  method = 'our_method';
%     write_auc_result_to_file([specie,'_',method],mic_auroc,mac_auroc,std_mac_auroc,specie);
