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
function [mic_auc,count,mic_auc_l] = evaluationYpredX(class_score, label)
if min(min(label))<0
    label = (label+1)/2;
end

nlabel = size(label,2);
nnode=size(label,2);
npos=zeros(1,nlabel);
st=[3,11,31,101,3];
ed=[10,30,100,300,300];
for i=1:nlabel
    npos(i)=nnz(label(:,i));
end
for j=1:5
    cat=npos(npos>=st(j) & npos<=ed(j));
    mic_auc = 0;
    count = 0;
    ap = 0;
    apc = 0;
    mic_auc_l=[];
    ap_l=[];
    for k=cat
        auc = calcROCarea(class_score(:,k),label(:,k));
        prec = calPR(class_score(:,k),label(:,k));
        if ~isnan(prec)
            ap = ap + prec;
            apc = apc + 1;
            ap_l=[ap_l,prec];
        end
        if ~isnan(auc)
            count = count + 1;
            mic_auc = mic_auc + auc;
            mic_auc_l=[mic_auc_l,auc];
        end
        class_score(:,k) = class_score(:,k);
        
    end
    s=class_score(:,cat);
    s=s(:);
    l=label(:,cat);
    l=l(:);
    mac_auc = calcROCarea(s,l);
    
    fprintf('%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\n',i,j,ap/apc,std(ap_l),mic_auc/count,std(mic_auc_l),mac_auc,length(cat));
    
end


end
