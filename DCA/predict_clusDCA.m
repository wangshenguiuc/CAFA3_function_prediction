function mac_auroc = predict_clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation,target_gene)
%CLUSDCA Summary of this function goes here
%   Detailed explanation goes here

Gene_GO_train_annotation = double(Gene_GO_train_annotation>0);

ly = GO_embedding;


active_gene = find(sum(Gene_GO_train_annotation,2)>0);
train_gene = setdiff(active_gene,target_gene);
test_gene = target_gene;
rng(2)

label = Gene_GO_train_annotation;
lx = Gene_embedding;
nnode = size(lx,1);
[nlabel,dim] = size(ly);
nfold = 3;
rp = randperm(nnode);
ntest = floor(nnode/nfold);

nrepeat = 10;

st_d=[1,11,31,101];
ed_d=[10,30,100,300];

final_score = zeros(nnode,nlabel);

npos =zeros(1,nlabel);
rp = randperm(nnode);
st = 1;
ed = ntest;
test_ind = test_gene;
train_ind = train_gene;
train_a = label;
train_a(test_ind,:) = 0;
label_tp = train_a;
w = embed(train_ind,label_tp,lx,ly);
class_score = lx*w*ly';
final_score(test_ind,:) = class_score(test_ind,:);
mac_auroc = mac_auc_evaluation(final_score(test_ind,:),label(test_ind,:));



end

