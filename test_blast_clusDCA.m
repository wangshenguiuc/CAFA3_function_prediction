addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
specie_l = {'Human','Yeast','Elegans','Drosophila','Mouse'};
specie_nick_name_l = {'Human','Yeasx','CAEEL','Drome','Rat'};
taxid_l = {'9606','4932','6239','7227','10090'};
for sp_ind =5:-1:4
    specie = specie_l{sp_ind};
    specie_nick_name = specie_nick_name_l{sp_ind};
    taxid = taxid_l{sp_ind};
    [GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
    [Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
    Gene_GO_train_annotation = read_annotation(specie,Gene_name,GO_name,GO_net);
    ngene = length(Gene_name);
    term_eia = pfp_eia(GO_net, logical(Gene_GO_train_annotation));
        
    ngene = length(Gene_name);
    
    GO_embedding_file_name = '../data/embedding_vector/clusDCA/GO_dim2500 rsp_0.8.US';
    GO_embedding = dlmread(GO_embedding_file_name);%learn_DCA_vector(GO_net,0.8,2500,GO_name,specie);
    embedding_file_name = ['../data/embedding_vector/clusDCA/',specie,'_Gene_dim2500 rsp_0.5.US'];
    Gene_embedding_base = dlmread(embedding_file_name);
% fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_MF_BP_CCloo_H.ant_M.ant_D.ant_E.ant.node.txt_500_10_0.500000_400_ppi_rwr_0_seq_rwr_0_enum_9_lr_0.025000_mode_2_meta__Y7E9G_Y7H9G_Y7M9G_Y7D9G_Y9G8G_Y9G';
% our_embedding_file_name = ['\\?\E:\swang141\project\SequencingNetwork\Sheng\data\embedding_vector\LINE\CV_ant\',fname];
%     Gene_embedding_our = read_embedding(Gene_name, our_embedding_file_name);
%     
    mode = 2;
    method ='Blast_clusDCA';
    afile = '..\data\annotation\blast_annotation\swiss_prot_go_annotation.txt';
    oa = read_blast_annotation(afile,GO_name,GO_net,mode,specie_nick_name);
    
    ifile = ['..\data\sequence_data\',taxid,'_to_uniprot.evalue.txt'];
    B = pfp_importblastp(ifile);
    
    qseqid = B.qseqid;
    pred = pfp_blast(qseqid(1:100),B,oa);    
    nfold = 3;    
    rng(2)
    [nnode,nlabel] = size(Gene_GO_train_annotation);
    ntest = floor(nnode/nfold);
    blast_final_score=zeros(nnode,nlabel);
    
    for p=1:nrepeat
        npos=zeros(1,nlabel);
        rp = randperm(nnode);
        st = 1;
        ed = ntest;
        test_ind = rp(st:ed);
        train_ind = rp(ed+1:nnode);
        
        test_oa = oa;
        
        [hold_set,qseqid,valid_pos] = find_hold_out_gene_name(values(Gene_name_rev,num2cell(test_ind)),oa.object,taxid);
        test_oa.object(hold_set) = [];
        test_oa.annotation(hold_set,:) = [];
        
        pred = pfp_blast(qseqid,B,oa);
        blast_final_score(test_ind(valid_pos),:) = pred.score;
        
       [~,~,~,~,clusDCA_final_score_base]=clusDCA(Gene_embedding_base, GO_embedding, Gene_GO_train_annotation, term_eia,GO_namespace,specie,GO_net,10);
         
%         [~,~,~,~,clusDCA_final_score_our]=clusDCA(Gene_embedding_our, GO_embedding, Gene_GO_train_annotation, term_eia,GO_namespace,specie,GO_net,10);
%         final_score_our = zscore(clusDCA_final_score_our)+ zscore(blast_final_score);
%         [mic_auroc,mac_auroc] = evaluation(final_score_our(test_ind,:),Gene_GO_train_annotation(test_ind,:),Gene_GO_train_annotation,GO_namespace,1);
        final_score_base = zscore(blast_final_score);
       [mic_auroc,mac_auroc] = evaluation(final_score_base(test_ind,:),Gene_GO_train_annotation(test_ind,:),Gene_GO_train_annotation,GO_namespace,1);
      
    end
    write_auc_result_to_file([specie,'_',method],avg_mic_auroc,avg_mac_auroc);
end