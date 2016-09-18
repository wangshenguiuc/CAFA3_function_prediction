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
function [mac_auroc] = mac_auc_evaluation(class_score, label)

nlabel = size(label,2);

st=[3,11,31,101,3];
ed=[10,30,100,300,300];
bias = 1;

for j=5:5
    cat=1:nlabel;
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
    mac_auroc= mac_auc/count;
end

end
