function [eval_res,avg_mic_auroc,avg_mac_auroc] = multiClusDCA_Blast(Gene_embedding, GO_embedding, Gene_GO_annotation, GO_namespace, test_sid,other_gene_go_annotation,Seq_net)
%CLUSDCA Summary of this function goes here
%   Detailed explanation goes here

lx = Gene_embedding{test_sid};
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
st_d = 1:10:50;
ed_d = 11:10:51;
% ed_d=[10,30,100,300];
eval_res = cell(1,nrepeat);
mic_auroc = cell(1,nrepeat);
mac_auroc = cell(1,nrepeat);
avg_mic_auroc = zeros(3,5);
avg_mac_auroc = zeros(3,5);
eval_res = [];
for p=1:nrepeat
    %     fprintf('start iteraction:%d / %d\n',p,nrepeat);
    npos=zeros(1,nlabel);
    rp = randperm(nnode{test_sid});
    for i=2:2
        func_t=GO_namespace(GO_namespace(:,2)==i,1);
        for tp=1:length(st_d)
            %         fprintf('case %d type %d / 3... \n',i, tp);
            st = 1;
            ed = ntest;
            test_ind = rp(st:ed);
            train_ind = rp(ed+1:nnode{test_sid});
            train_a = Gene_GO_annotation{test_sid};
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
                ct = ct+1;
                for g=func_tp'
                    blast_score(tg, g) = max(Seq_net(tg,go_anot_gene{g}));
                end
            end            
        end
    end
    
    final_score = blast_score + embed_score;
    
    [mic_auroc{p},mac_auroc{p}] = evaluation(final_score(test_ind,:),Gene_GO_annotation{test_sid}(test_ind,:),Gene_GO_annotation{test_sid},GO_namespace,1);
    
    avg_mic_auroc = avg_mic_auroc + mic_auroc{p}/nrepeat;
    avg_mac_auroc = avg_mac_auroc + mac_auroc{p}/nrepeat;
    
end

end

