function mac_auroc = cv_clusDCA(Gene_embedding, GO_embedding, Gene_GO_train_annotation)
%CLUSDCA Summary of this function goes here
%   Detailed explanation goes here

Gene_GO_train_annotation = double(Gene_GO_train_annotation>0);

ly = GO_embedding;


active_gene = find(sum(Gene_GO_train_annotation,2)>0);
rng(2)

label = Gene_GO_train_annotation(active_gene,:);
lx = Gene_embedding(active_gene,:);
nnode = size(lx,1);
[nlabel,dim] = size(ly);
nfold = 3;
rp = randperm(nnode);
ntest = floor(nnode/nfold);

nrepeat = 10;

st_d=[1,11,31,101];
ed_d=[10,30,100,300];
mac_auroc = zeros(nrepeat,1);

final_score = zeros(nnode,nlabel);
for p=1:nrepeat
    fprintf('start iteration:%d / %d\n',p,nrepeat);
    npos=zeros(1,nlabel);
    rp = randperm(nnode);
    st = 1;
    ed = ntest;
    test_ind = rp(st:ed);
    train_ind = rp(ed+1:nnode);
    train_a = label;
    train_a(test_ind,:) = 0;
    label_tp = train_a;
    w = embed(train_ind,label_tp,lx,ly);
    %     w = nonlinear_sigmodSgd(a,label,test_ind,train_ind,lx,ly,D,w);
    class_score = lx*w*ly';
    final_score(test_ind,:) = class_score(test_ind,:);
    mac_auroc(p,:) = mac_auc_evaluation(final_score(test_ind,:),label(test_ind,:));
    fprintf('%d %f',p,mac_auroc(p,:));
end

mac_auroc = mean(mac_auroc);

end

