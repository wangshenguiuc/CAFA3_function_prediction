function [eval_res,avg_mic_auroc,avg_mac_auroc] = multiClusDCA(Gene_embedding, GO_embedding, Gene_GO_annotation, term_eia,GO_namespace, test_sid,GO_net,SortDist)
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
final_score=zeros(nnode{test_sid},nlabel);

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
            multi_label_tp = train_a(:,func_tp);
            multi_lx = lx;
            weight = ones(1,nnode{test_sid});
            for sp=[2]
                if sp == test_sid
                    continue;
                end
                select_ind = SortDist{sp}.I(test_ind,1:1);
                select_ind = unique(select_ind(:));
                multi_lx = [multi_lx;Gene_embedding{sp}(select_ind,:)];
                multi_label_tp = [multi_label_tp;Gene_GO_annotation{sp}(select_ind,func_tp)];
                weight = [weight,0.5*ones(1,length(select_ind))];
%                 length(select_ind)
            end
            
            all_node = size(multi_lx,1);
            multi_ind = setdiff(1:all_node,test_ind);
            multi_ly = ly(func_tp,:);
            w1 = embed(multi_ind,multi_label_tp,multi_lx,multi_ly,weight);
            w2 = embed(train_ind,train_a(:,func_tp),lx,multi_ly,ones(1,nnode{test_sid}));
            class_score = lx*w1*ly(func_tp,:)';
            final_score(test_ind,func_tp) = class_score(test_ind,:);
        end
    end
    %     final_score(test_ind,:) = optimize_score_GO_propogation(final_score(test_ind,:),GO_net);
    %     for i=1:nlabel
    %         final_score(:,i) = final_score(:,i) - mean(final_score(:,i));
    %     end
    %     fprintf('final evaluation\n');
    [mic_auroc{p},mac_auroc{p}] = evaluation(final_score(test_ind,:),Gene_GO_annotation{test_sid}(test_ind,:),Gene_GO_annotation{test_sid},GO_namespace,1);
    %          [mic_auroc{p},mac_auroc{p}] = evaluation(final_score(test_ind,:),Gene_GO_annotation(test_ind,:),Gene_GO_annotation,GO_namespace,0)
    
    %     eval_res{p}.term_eval_res = term_evaluation(final_score(test_ind,:),Gene_GO_annotation(test_ind,:),Gene_GO_annotation,GO_namespace,1);
    
    %    eval_res{p}.protein_eval_res = protein_evaluation(final_score(test_ind,:),Gene_GO_annotation(test_ind,:),term_eia,Gene_GO_annotation,GO_namespace);
    avg_mic_auroc = avg_mic_auroc + mic_auroc{p}/nrepeat;
    avg_mac_auroc = avg_mac_auroc + mac_auroc{p}/nrepeat;
    
end

end

