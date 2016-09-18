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
function [mic_auroc,mac_auroc,mac_auroc_detail,std_mac_auroc] = evaluation(class_score, label,global_label,GO_namespace,bias)

g1 = GO_namespace(:,1);
g2 = GO_namespace(:,2);
nlabel = size(label,2);
go_type = unique(g2);
npos=zeros(1,nlabel);
st=[3,11,31,101,3];
ed=[10,30,100,300,300];
for i=1:nlabel
    npos(i)=nnz(global_label(:,i));
end
mic_auroc = zeros(2,4);
mac_auroc = zeros(2,4);
std_mac_auroc = zeros(2,4);
mac_auroc_detail = cell(2,4);
for i=[1,2]%go_type'
    func_t=g1(g2==i);
    for j=1:4
        cat=func_t(npos(func_t)>=st(j) & npos(func_t)<=ed(j));
        cat = cat';
        count = 0;
        mac_auc_l=[];
        mac_auc=0;
        for k=cat
            auc = calcROCarea(class_score(:,k),label(:,k));
            if ~isnan(auc)
                count = count + 1;
                mac_auc = mac_auc + auc;
                mac_auc_l=[mac_auc_l,auc];
            end
            if bias==1
                class_score(:,k) = class_score(:,k)-mean(class_score(:,k));
            end
        end
        s=class_score(:,cat);
        s=s(:);
        l=label(:,cat);
        l=l(:);
        mic_auc = calcROCarea(s,l);
        mic_auroc(i,j) = mic_auc;
        mac_auroc(i,j) = mac_auc/count;
        std_mac_auroc(i,j) = std(mac_auc_l);
        mac_auroc_detail{i,j}  = mac_auc_l;
        fprintf('type=%d\tcat=%d\tmacro AUC=%f\tstd(macro AUC)=%f\tmicro AUC=%f\tcount=%d\n',i,j,mac_auc/count,std(mac_auc_l),mic_auc,length(cat));
        %
    end
    
end
mic_auroc = mic_auroc';
mic_auroc = mic_auroc(:);
mac_auroc = mac_auroc';
mac_auroc = mac_auroc(:);
std_mac_auroc = std_mac_auroc';
std_mac_auroc = std_mac_auroc(:);
end
