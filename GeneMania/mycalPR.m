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
function [ap] = calPR(score, label)

nnpos=nnz(label);
[o,v]=sort(score,'descend');
cor=0;
tot = 0;
ap = 0;
for i=1:length(v)
    id = v(i);
    tot = tot + 1;
    if label(id)==1
        cor = cor +1;
        ap = ap + cor/tot;
    end 
    if cor==ceil(0.1*nnpos)
        break
    end
end
ap = ap/cor;
end
