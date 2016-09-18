function [Gene_GO_train_annotation,tgt_gene_id,Gene_embedding,Gene_name,Gene_name_rev]  = ...
read_CAFA2_annotation(Gene_name,GO_name,GO_net,Gene_embedding,Gene_name_rev)
addpath '../data/annotation/CAFA2'
 addpath '../data/CAFA2_target_gene'

[e1,e2,~,~] = textread('hpoa.annotation.t0','%s%s%s%s');

filt = GO_name.isKey(e1) & Gene_name.isKey(e2);
e1 = e1(filt);
go = cell2mat(values(GO_name,e1));

e2 = e2(filt);
gene = cell2mat(values(Gene_name,e2));

GO_Gene_mat = sparse(go,gene,1,size(GO_name,1),size(Gene_name,1));

old_GO_Gene_mat = GO_Gene_mat;

it = 1;
while it==1 || ~isequal(old_GO_Gene_mat , GO_Gene_mat)
    old_GO_Gene_mat = GO_Gene_mat;
    GO_Gene_mat = GO_net * GO_Gene_mat + old_GO_Gene_mat;
%     fprintf('it:%d, edge num:%d\n',it,length(find(GO_Gene_mat(:)>0)));
    GO_Gene_mat = double(GO_Gene_mat>=1);
    it = it+1;
end

Gene_GO_train_annotation = GO_Gene_mat';

Gene_GO_train_annotation = double(Gene_GO_train_annotation>0);


[e1,e2,~,~] = textread('hpoa.annotation.t1','%s%s%s%s');

filt = GO_name.isKey(e1) & Gene_name.isKey(e2);
e1 = e1(filt);
go = cell2mat(values(GO_name,e1));

e2 = e2(filt);
gene = cell2mat(values(Gene_name,e2));

GO_Gene_mat = sparse(go,gene,1,size(GO_name,1),size(Gene_name,1));

old_GO_Gene_mat = GO_Gene_mat;

it = 1;
while it==1 || ~isequal(old_GO_Gene_mat , GO_Gene_mat)
    old_GO_Gene_mat = GO_Gene_mat;
    GO_Gene_mat = GO_net * GO_Gene_mat + old_GO_Gene_mat;
%     fprintf('it:%d, edge num:%d\n',it,length(find(GO_Gene_mat(:)>0)));
    GO_Gene_mat = double(GO_Gene_mat>=1);
    it = it+1;
end

Gene_GO_train_annotation = GO_Gene_mat';

Gene_GO_train_annotation = double(Gene_GO_train_annotation>0);

% tgt_gene = textread('xxo_all_typex.txt','%s');
% filter = isKey(Gene_name,tgt_gene);
% avg_embedding = mean(Gene_embedding);
% ngene = length(Gene_name);
% impute_gene = find(filter==0);
% ct = 1;
% for i = impute_gene'
%     Gene_name_rev(ngene+ct) = char(tgt_gene(i));
%     Gene_name(char(tgt_gene(i))) = ngene+ct;
%     Gene_embedding = [Gene_embedding;avg_embedding];
%     ct = ct+1;
% end
% tgt_gene_id = cell2mat(values(Gene_name,tgt_gene));

