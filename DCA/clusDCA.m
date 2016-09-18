function [eval_res,mic_auroc,mac_auroc,mac_auroc_detail,final_score] = clusDCA(Gene_embedding, GO_embedding, Gene_GO_annotation, term_eia,GO_namespace, specie,GO_net,step_size)
%CLUSDCA Summary of this function goes here
%   Detailed explanation goes here

%
% specie='Human';
% dim_num=2500;
%     lx = dlmread(['/srv/data/swang141/swang141/OutputMatrix/NoIsoString',specie,'USX d=',num2str(dim),' us=1 rsp=0.5.txt']);
%     ly = dlmread(['/srv/data/swang141/swang141/OutputMatrix/NoIsoString',specie,'USY d=',num2str(dim),' us=1 rsp=0.8 bp=0.8.txt']);
lx = Gene_embedding;
ly = GO_embedding;
[nnode, Gene_dim] = size(Gene_embedding);
[nlabel,GO_dim] = size(GO_embedding);

nfold = 3;

rng(2)
% D = squareform(pdist(lx,'cosine'));

ntest = floor(nnode/nfold);

final_score=zeros(nnode,nlabel);

nrepeat = 1;
if step_size == -1
st_d=[1,11,31,101];
ed_d=[10,30,100,300];
else
st_d = 0:step_size:500;
ed_d = step_size:step_size:500+step_size;
end
eval_res = cell(1,nrepeat);
mic_auroc = zeros(8,nrepeat);
mac_auroc = zeros(8,nrepeat);
avg_mic_auroc = zeros(2,4);
avg_mac_auroc = zeros(2,4);
eval_res = [];
for p=1:nrepeat
    %     fprintf('start iteraction:%d / %d\n',p,nrepeat);
    npos=zeros(1,nlabel);
    rp = randperm(nnode);
    st = 1;
    ed = ntest;
    test_ind = rp(st:ed);
    train_ind = rp(ed+1:nnode);
    train_a = Gene_GO_annotation;
    train_a(test_ind,:) = 0;
    for k=1:nlabel
        npos(k)=nnz(train_a(:,k));
    end
    for i=1:2
        func_t=GO_namespace(GO_namespace(:,2)==i,1);
        
        for tp=1:length(st_d)            
            func_tp = func_t(npos(func_t)>=st_d(tp)*2/3& npos(func_t)<=ed_d(tp)*2/3);
            label_tp = train_a(:,func_tp);
            
            w = embed(train_ind,label_tp,lx,ly(func_tp,:),ones(1,nnode));
            class_score = lx*w*ly(func_tp,:)';
            final_score(test_ind,func_tp) = class_score(test_ind,:);
        end
    end
    %     final_score(test_ind,:) = optimize_score_GO_propogation(final_score(test_ind,:),GO_net);
    %     for i=1:nlabel
    %         final_score(:,i) = final_score(:,i) - mean(final_score(:,i));
    %     end
    %     fprintf('final evaluation\n');
    [mic_auroc(:,p),mac_auroc(:,p),mac_auroc_detail] = evaluation(final_score(test_ind,:),Gene_GO_annotation(test_ind,:),Gene_GO_annotation,GO_namespace,1);
    %          [mic_auroc{p},mac_auroc{p}] = evaluation(final_score(test_ind,:),Gene_GO_annotation(test_ind,:),Gene_GO_annotation,GO_namespace,0)
    
    %     eval_res{p}.term_eval_res = term_evaluation(final_score(test_ind,:),Gene_GO_annotation(test_ind,:),Gene_GO_annotation,GO_namespace,1);
    
    %    eval_res{p}.protein_eval_res = protein_evaluation(final_score(test_ind,:),Gene_GO_annotation(test_ind,:),term_eia,Gene_GO_annotation,GO_namespace);
%     avg_mic_auroc = avg_mic_auroc + mic_auroc{p}/nrepeat;
%     avg_mac_auroc = avg_mac_auroc + mac_auroc{p}/nrepeat;
    
end
% avg_mic_auroc = avg_mic_auroc';
% avg_mic_auroc = avg_mic_auroc(:);
% avg_mac_auroc = avg_mac_auroc';
% avg_mac_auroc = avg_mac_auroc(:);

end

