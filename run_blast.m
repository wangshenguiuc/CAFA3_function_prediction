addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
specie = 'Human';%'Mouse';%'Yeast';%'Elegans'; % CAEEL
specie_nick_name = 'Human';%'Rat';%'CAEEL';%'Drome';
taxid = '9606';%'10090';%'6239'; %7227
[GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
[Gene_name,Gene_name_rev,Gene_net] = read_string_network(specie);
Gene_GO_train_annotation = read_annotation(specie,Gene_name,GO_name,GO_net);
ngene = length(Gene_name);
term_eia = pfp_eia(GO_net, logical(Gene_GO_train_annotation));
method_l={'Blast_all','Blast_1_specie','Blast_5_species'};
for mode = [0]
    method = method_l{mode+1};
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
    final_score=zeros(nnode,nlabel);
    nrepeat = 3;
    st_d=[1,11,31,101];
    ed_d=[10,30,100,300];
    eval_res = cell(1,nrepeat);
    mic_auroc = cell(1,nrepeat);
    mac_auroc = cell(1,nrepeat);
    avg_mic_auroc = [];
avg_mac_auroc = [];
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
        final_score(test_ind(valid_pos),:) = pred.score;
        [mic_auroc,mac_auroc] = evaluation(final_score(test_ind,:),Gene_GO_train_annotation(test_ind,:),Gene_GO_train_annotation,GO_namespace,0);
        
    avg_mic_auroc = [avg_mic_auroc, mic_auroc];
avg_mac_auroc = [avg_mac_auroc, mac_auroc];
    end
    write_auc_result_to_file([specie,'_',method],avg_mic_auroc,avg_mac_auroc);
    
end

