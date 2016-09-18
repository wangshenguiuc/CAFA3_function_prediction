function final_score = test_clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation,Gene_GO_test_annotation, GO_namespace,specie,human_target_id)
%CLUSDCA Summary of this function goes here
%   Detailed explanation goes here
rng(2)

ntgt = length(human_target_id);
ngene = size(Gene_GO_train_annotation,1);
ngo = size(Gene_GO_train_annotation,2);

active_gene = 1:ngene;%setdiff(1:ngene,human_target_id);
ly = GO_embedding;

label = Gene_GO_train_annotation(active_gene,:);

lx = Gene_embedding(active_gene,:);
nnode = size(lx,1);

pos_ind=zeros(ngo,ngene);
pos_l = zeros(ngo,1);
neg_ind=zeros(ngo,ngene);
neg_l = zeros(ngo,1);

for i=1:ngo
    p = find(Gene_GO_train_annotation(:,i)==1);
    n = find(Gene_GO_train_annotation(:,i)==0);
    pos_ind(i,1:length(p)) = p;
    pos_l(i) = length(p);
    neg_ind(i,1:length(n)) = n;
    neg_l(i) = length(n);
end

final_score=zeros(ntgt,ngo);

st_d=[1,11,31,101];
ed_d=[10,30,100,300];
func_t = 1:ngo;
for tp=1:4
    
    train_a = Gene_GO_train_annotation;
    for k=1:ngo
        npos(k)=nnz(train_a(:,k));
    end
    func_tp = func_t(npos(func_t)>=st_d(tp) & npos(func_t)<=ed_d(tp));
    label_tp = train_a(:,func_tp);
    
    w = embed(1:nnode,label_tp,lx,ly(func_tp,:));
    class_score = lx*w*ly(func_tp,:)';
    final_score(:,func_tp) = class_score(human_target_id,:);
end

[mic_auroc,mac_auroc] = evaluation(final_score,Gene_GO_test_annotation(human_target_id,:),Gene_GO_test_annotation,GO_namespace)

end

