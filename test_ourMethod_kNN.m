% addpath 'DCA/'
% addpath 'util/'
% addpath 'readData/'
% addpath 'LINE/'
% addpath 'blast/'
% specie_l = {'Human','Yeast','Drosophila','Mouse','Elegans'};
% test_sid = 1;
% specie = specie_l{test_sid};
% [GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
% Gene_name = cell(1,5);
% Gene_name_rev = cell(1,5);
% Gene_net = cell(1,5);
% Gene_GO_train_annotation = cell(1,5);
% Gene_embedding = cell(1,5);
% 
% GO_embedding_file_name = '../data/embedding_vector/clusDCA/GO_dim2500 rsp_0.8.US';
% GO_embedding = dlmread(GO_embedding_file_name);
% 
% fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_GOloo_H_M.asc_Y.asc_D.asc_E.asc.node.txt_500_10_0.500000_50_ppi_rwr_1_seq_rwr_1_enum_9_lr_0.025000_mode_1_meta__H7M9G_H7E9G_H7Y9G_H7D9G_H9G8G_H9G';
% embedding_file_name =['\\?\E:\swang141\project\SequencingNetwork\Sheng\data\embedding_vector\selected_code\',specie,'\',fname];
% 
% for sp=1:5
%     spt = char(specie_l(sp));
%     [Gene_name{sp},Gene_name_rev{sp}] = read_string_network(spt);
%     Gene_GO_train_annotation{sp} = read_annotation(spt,Gene_name{sp},GO_name,GO_net);
%     Gene_embedding{sp} = read_embedding(Gene_name{sp}, embedding_file_name);
% end
% 
% 
% other_gene_go_annotation = [];
% other_embedding = [];
% for sp=1:5
%     if sp==test_sid
%         continue
%     end
%     other_gene_go_annotation = [other_gene_go_annotation;Gene_GO_train_annotation{sp}];
%     other_embedding = [other_embedding;Gene_embedding{sp}];
% end
% 
% lx = Gene_embedding{test_sid};
% nspecie = length(Gene_embedding);
% nlabel = size(GO_embedding,1);
% nnode = cell(1,nspecie);
% for i=1:nspecie
%     nnode{i}= size(Gene_embedding{i},1);
% end
% 
% nfold = 3;
% 
% rng(2)
% 
% ntest = floor(nnode{test_sid}/nfold);
% 
% 
% nrepeat = 1;
% st_d = 1:10:501;
% ed_d = 11:10:511;
% npos=zeros(1,nlabel);
% rp = randperm(nnode{test_sid});
% st = 1;
% ed = ntest;
% test_ind = rp(st:ed);
% train_ind = rp(ed+1:nnode{test_sid});
% 
% test_embedding = lx(test_ind,:);
% same_spe_embedding = lx;
% other_spe_embedding = other_embedding;
% fprintf('start to get distnace matrix');
% same_D = pdist2(test_embedding,same_spe_embedding, 'cosine');
% other_D = pdist2(test_embedding,other_spe_embedding, 'cosine');
% same_ant = Gene_GO_train_annotation{test_sid};
% same_ant(test_ind,:) = 0;
other_score=zeros(nnode{test_sid},nlabel);
final_score=zeros(nnode{test_sid},nlabel);
same_score=zeros(nnode{test_sid},nlabel);
% other_ant = other_gene_go_annotation;
for nvote=1000:2000:11000
for cat=1:2
    func_t=GO_namespace(GO_namespace(:,2)==cat,1);
%         for i = 1:ntest
%             [v, o] = sort(same_D(i, train_ind));
%             o = o(~isinf(v));
%             k = min(nvote, length(o));
%             votes = sum(bsxfun(@rdivide, same_ant(train_ind(o(1:k)),func_t)', same_D(i, train_ind(o(1:k)))), 2);
%             same_score(test_ind(i),func_t) = votes;
%         end
        for i = 1:ntest
            [v, o] = sort(other_D(i, train_ind));
            o = o(~isinf(v));
            k = min(nvote, length(o));
            votes = sum(bsxfun(@rdivide, other_ant(train_ind(o(1:k)),func_t)', other_D(i, train_ind(o(1:k)))), 2);
            other_score(test_ind(i),func_t) = votes;
        end    
end
nvote
        final_score =other_score;
fprintf('start other\n')
[mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc]  = evaluation(final_score(test_ind,:),Gene_GO_train_annotation{test_sid}(test_ind,:),Gene_GO_train_annotation{test_sid},GO_namespace,1);

% final_score = same_score;
% fprintf('start same\n')
% [mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc]  = evaluation(final_score(test_ind,:),Gene_GO_train_annotation{test_sid}(test_ind,:),Gene_GO_train_annotation{test_sid},GO_namespace,1);

final_score = other_score+same_score;
fprintf('start both\n')
[mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc]  = evaluation(final_score(test_ind,:),Gene_GO_train_annotation{test_sid}(test_ind,:),Gene_GO_train_annotation{test_sid},GO_namespace,1);
end