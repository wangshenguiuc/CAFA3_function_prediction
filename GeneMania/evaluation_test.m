% Evaluate prediction performance.
%
% [Input]
% class_score: (# test case) x (# class) matrix of predicted class scores.
%              Higher number represents higher confidence.
% label: (# test case) x (# class) matrix of ground truth annotations.
%        Note each case can have multiple labels
%
% [Output]
% acc: accuracy measure. Pick top prediction for each test case and
%      see how often it matches with one of the true labels.
% f1: micro-averaged F1-score. Pick top alpha predictions for each test case,
%     calculate the contigency table for each class, sum up the table across
%     all classes then calculate the F1-score.
%
function [mic_auc,mac_auc,ap_10,ap] = evaluation_test(class_score, label,global_label,specie,i,j,label_func,remove_list)
if min(min(label))<0
    label = (label+1)/2;
end
if min(min(global_label))<0
    global_label = (global_label+1)/2;
end
if strcmp(specie,'Mouse')
    selectGO = textread(['/home/swang141/research/Bio Network Ontology Prediction/GoPrediction/Data/',specie,'Graph/selectgo.txt'], '%d');
end
[g1,g2] = textread(['/home/swang141/research/Bio Network Ontology Prediction/GoPrediction/Data/',specie,'Graph/noisogo_namespace.txt'], '%d%d');
nlabel = size(label,2);
npos=zeros(1,nlabel);
st=[3,11,31,101,3];
ed=[10,30,100,300,300];
for ii=1:nlabel
    npos(ii)=nnz(global_label(:,ii));
end
% npos(1:5)
mac_auc_l=[];
% for i=1:nnode
%     class_score(i,:) = class_score(i,:) - max(class_score(i,:));
% end
func_t=g1(g2==i);

cat=func_t(npos(func_t)>=st(j) & npos(func_t)<=ed(j));
cat = label_func;
cat = setdiff(cat,remove_list);
if strcmp(specie,'Mouse')
    cat = intersect(cat,selectGO);
end
mac_auc = 0;
count = 0;
ap = 0;
apc = 0;
apc_10=0;
ap_10=0;
mac_auc_l=[];
ap_l=[];
for k=cat
    
    auc = calcROCarea(class_score(:,k),label(:,k));
% %                 ap_at_10=0;
%     prec = 0;
    
    if nnz(label(:,k))>0
        ap_at_10 = calcPR(class_score(:,k),label(:,k));
        ap_10 = ap_10 + ap_at_10;
        apc_10 = apc_10 + 1;
    end
    if nnz(label(:,k))>0
        [Xpr,Ypr,Tpr,prec] = perfcurve(full(label(:,k)), full(class_score(:,k)), 1, 'xCrit', 'reca', 'yCrit', 'prec');
        ap = ap + prec;
        apc = apc + 1;
    end
    if ~isnan(auc)
        count = count + 1;
        mac_auc = mac_auc + auc;
    end
    %             class_score(:,k) = class_score(:,k)-mean(class_score(:,k));
    
end
s=class_score(:,cat);
s=s(:);
l=label(:,cat);
l=l(:);
mic_auc = calcROCarea(s,l);
% micprec = 0;
% [~,~,~,micprec] = perfcurve(l, s, 1, 'xCrit', 'reca', 'yCrit', 'prec');
% fprintf('%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n',i,j,ap_10/apc_10,ap/apc,std(ap_l),mic_auc/count,std(mic_auc_l),mac_auc,length(cat));
% count
% length(cat)
% length(remove_list)
mac_auc = mac_auc/count;
ap_10 = ap_10/apc_10;
ap = ap/apc;
end
