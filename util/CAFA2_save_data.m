function [pred] = CAFA2_save_data( final_score, tp, CAFA_gene_id)
%CAFA2_SAVE_DATA Summary of this function goes here
%   Detailed explanation goes here
 addpath '../data/CAFA2_target_gene'
load(['..\analysis\CAFA2\github\CAFA2-master\baselines\',tp,'\BB4H.mat']);
tgt_id = textread('xxo_all_typex.txt','%s');


pred.score = zeros(size(pred.score ));
pred.score(CAFA_gene_id,:) = final_score;

pred.score = pred.score-min(pred.score(:));

pred.score = pred.score/max(pred.score(:));
pred_score = pred.score;
save(['..\analysis\CAFA2\github\CAFA2-master\baselines\',tp,'\our_most_popular'],'pred');

% tgt_id = textread('xxo_all_typex.txt','%s');
% 
% ind = [];
% for t = pred.object'
%     ind = [ind,find(strcmp(t,tgt_id))];    
% end
% 
% pred.score = final_score(ind,:);
% 
% pred.score = pred.score-min(pred.score(:));
% 
% pred.score = pred.score/max(pred.score(:));
% % pred.score = repmat(final_score,1350, 1);
% save(['..\analysis\CAFA2\github\CAFA2-master\baselines\',tp,'o\our_most_popular'],'pred');

% load(['..\analysis\CAFA2\github\CAFA2-master\benchmark\groundtruth\',tp,'oa.mat']);

% mac_auroc=mac_auc_evaluation(final_score(tgt_edge_ct>10,:), oa.annotation(tgt_edge_ct>0,:))