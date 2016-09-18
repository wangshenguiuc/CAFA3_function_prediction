function pr= calcPR_v2(scores, y,recallPoint)
%function PR50 = calcPR(Scores, Labels)
% calculates precision at a given recall level
%  Datapoints are the rows, labelling categories are the columns
%  y is predicted scores, t is label vector
%  Positives are labelled with ones, negatives are labelled with zeros.

pos = find(y==1);

for ii = 1:length(pos)
    thresh(ii) = scores(pos(ii));
    recall(ii) = length(intersect(pos,find(scores>=thresh(ii))));
    precision(ii) = length(intersect(pos,find(scores>=thresh(ii))))./length(find(scores>=thresh(ii)));
end
TP = sum(y==1);
if(TP<1);
    pr = NaN;
else

recall = recall./TP;
diff = recall-recallPoint;
ss = find(diff<0);
diff(ss) = inf;
mymin = find(diff == min(diff));
pr = precision(mymin);
pr = min(pr);


end