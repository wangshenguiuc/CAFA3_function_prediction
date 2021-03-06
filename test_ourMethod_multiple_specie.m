addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
specie_l = {'Human','Yeast','Drosophila','Mouse','Elegans'};
test_sid = 1;

[GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
Gene_name = cell(1,5);
Gene_name_rev = cell(1,5);
Gene_net = cell(1,5);
Gene_GO_train_annotation = cell(1,5);
Gene_embedding = cell(1,5);

GO_embedding_file_name = '../data/embedding_vector/clusDCA/GO_dim2500 rsp_0.8.US';
GO_embedding = dlmread(GO_embedding_file_name);

for sp=1:5
    specie = char(specie_l(sp));
    [Gene_name{sp},Gene_name_rev{sp}] = read_string_network(specie);
    Gene_GO_train_annotation{sp} = read_annotation(specie,Gene_name{sp},GO_name,GO_net);
    %     Gene_embedding{sp} = read_embedding(Gene_name{sp}, embedding_file_name);
end
specie = specie_l{test_sid};
fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_GOloo_H_M.asc_Y.asc_D.asc_E.asc.node.txt_500_10_0.500000_50_ppi_rwr_1_seq_rwr_1_enum_9_lr_0.025000_mode_1_meta__H7M9G_H7E9G_H7Y9G_H7D9G_H9G8G_H9G';
embedding_file_name =['\\?\E:\swang141\project\SequencingNetwork\Sheng\data\embedding_vector\selected_code\',specie,'\',fname];
Our_Gene_embedding = read_embedding(Gene_name{test_sid}, embedding_file_name);

other_gene_go_annotation = [];
Seq_net = [];
for sp=1:5
    if sp==test_sid
        continue
    end
    other_gene_go_annotation = [other_gene_go_annotation;Gene_GO_train_annotation{sp}];
    %     Seq_net = [Seq_net,Gene_embedding{test_sid}*Gene_embedding{sp}'];
end

Seq_net = sparse(zeros(length(Gene_name{test_sid}),size(other_gene_go_annotation,1)));
[Seq_net_mul,Seq_Gene_name,Seq_Gene_name_rev] = read_sequence_network();

test_gname = keys(Gene_name{test_sid});
filter1 = isKey(Seq_Gene_name,test_gname);
test_gname = test_gname(filter1);
test_gname_seq_id = cell2mat(values(Seq_Gene_name,test_gname));

other_gname_name = [];
for sp=1:5
    if sp==test_sid
        continue
    end
    other_gname_name = [other_gname_name,keys(Gene_name{sp})];
end
filter2 = isKey(Seq_Gene_name,other_gname_name);
other_gname_name = other_gname_name(filter2);
other_gname_seq_id = cell2mat(values(Seq_Gene_name,other_gname_name));

Seq_net_mul_new = sparse(zeros(size(Seq_net)));
Seq_net_mul_new(filter1,filter2) = Seq_net_mul(test_gname_seq_id,other_gname_seq_id);
clear Seq_net_mul;
clear Seq_net;
embedding_file_name = ['../data/embedding_vector/clusDCA/',specie_l{test_sid},'_Gene_dim2500 rsp_0.5.US'];
Gene_embedding{test_sid} = dlmread(embedding_file_name);

% lx = Gene_embedding{test_sid};
lx = Our_Gene_embedding;
ly = GO_embedding;
nspecie = length(Gene_embedding);
nlabel = size(GO_embedding,1);
nnode = cell(1,nspecie);
for i=1:nspecie
    nnode{i}= size(Gene_embedding{i},1);
end

nfold = 3;

rng(2)
% D = squareform(pdist(lx,'cosine'));

ntest = floor(nnode{test_sid}/nfold);
embed_score=zeros(nnode{test_sid},nlabel);
final_score=zeros(nnode{test_sid},nlabel);
blast_score=zeros(nnode{test_sid},nlabel);

nrepeat = 1;
% st_d=[1,11,31,101];
st_d = 1:10:501;
ed_d = 11:10:511;
% ed_d=[10,30,100,300];
eval_res = cell(1,nrepeat);
mic_auroc = cell(1,nrepeat);
mac_auroc = cell(1,nrepeat);
avg_mic_auroc = zeros(3,5);
avg_mac_auroc = zeros(3,5);
eval_res = [];
%     fprintf('start iteraction:%d / %d\n',p,nrepeat);
npos=zeros(1,nlabel);
rp = randperm(nnode{test_sid});
for i=1:2
    func_t=GO_namespace(GO_namespace(:,2)==i,1);
    for tp=1:length(st_d)
        %         fprintf('case %d type %d / 3... \n',i, tp);
        st = 1;
        ed = ntest;
        test_ind = rp(st:ed);
        train_ind = rp(ed+1:nnode{test_sid});
        train_a = Gene_GO_train_annotation{test_sid};
        train_a(test_ind,:) = 0;
        for k=1:nlabel
            npos(k)=nnz(train_a(:,k));
        end
        func_tp = func_t(npos(func_t)>=st_d(tp)*2/3& npos(func_t)<=ed_d(tp)*2/3);
        multi_ly = ly(func_tp,:);
        w = embed(train_ind,train_a(:,func_tp),lx,multi_ly,ones(1,nnode{test_sid}));
        class_score = lx*w*ly(func_tp,:)';
        embed_score(test_ind,func_tp) = class_score(test_ind,:);
        
        
        go_anot_gene = cell(1,nlabel);
        active_label = [];
        for g=func_tp'
            go_anot_gene{g} = find(other_gene_go_annotation(:,g)>0);
        end
        for tg=test_ind
            for g=func_tp'
                if isempty(go_anot_gene{g})
                    continue;
                end
                blast_score(tg, g) = max(Seq_net_mul_new(tg,go_anot_gene{g}));
            end
        end
    end
end
blast_score_norm = (blast_score-min(blast_score(:)))/(max(blast_score(:))-min(blast_score(:)));
embed_score_norm = (embed_score-min(embed_score(:)))/(max(embed_score(:))-min(embed_score(:)));
final_score =  zscore(blast_score')' +zscore(embed_score_norm')';
[mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc]  = evaluation(final_score(test_ind,:),Gene_GO_train_annotation{test_sid}(test_ind,:),Gene_GO_train_annotation{test_sid},GO_namespace,1);
method = 'Blast_clusDCA_zs';
% write_auc_result_to_file([specie,'_',method],mic_auroc,mac_auroc,std_mac_auroc);

final_score =  blast_score_norm+embed_score_norm;
[mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc]  = evaluation(final_score(test_ind,:),Gene_GO_train_annotation{test_sid}(test_ind,:),Gene_GO_train_annotation{test_sid},GO_namespace,1);
method = 'Blast_clusDCA_ratio';
% write_auc_result_to_file([specie,'_',method],mic_auroc,mac_auroc,std_mac_auroc);

final_score =  blast_score_norm ;
[mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc]  = evaluation(final_score(test_ind,:),Gene_GO_train_annotation{test_sid}(test_ind,:),Gene_GO_train_annotation{test_sid},GO_namespace,1);
method = 'Blast_5specie_select_code';
% write_auc_result_to_file([specie,'_',method],mic_auroc,mac_auroc,std_mac_auroc);


