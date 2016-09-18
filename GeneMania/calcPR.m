function [PR50,FP] = calcPR(y, t,recallPoint)
%function PR50 = calcPR(Scores, Labels)
% calculates precision at a given recall level
%  Datapoints are the rows, labelling categories are the columns
%  y is predicted scores, t is label vector
%  Positives are labelled with ones, negatives are labelled with zeros.

t = t > 0 ;
y = full(y);
[ysorted,indx]=sort(y,'descend');
[ysortedR,idxR] = sort(y(end:-1:1),'descend');
idxR = length(y)-idxR + 1;
myt = t;
t = t(indx);
t2 = myt(idxR);

numElem = 1:length(t);

recall = cumsum(t)/sum(t);
recall2 = cumsum(t2)/sum(t2);


my_tp = cumsum(t);
my_fp = cumsum(~t);

fp = cumsum(~t)/sum(~t);

perc = cumsum(t)./(cumsum(~t)+cumsum(t));

perc2 = cumsum(t2)./(cumsum(~t2)+cumsum(t2));


[junk,ix] = (min(abs(recall-recallPoint)));
[junk2,ix2] = min(abs(recall2-recallPoint));
%ix = ix(end);
PR1 = (perc(ix(1))); 
PR2 = perc2(ix2(1));

r2 = recall-recallPoint;
r2(r2<0) = inf;
[junk,ix] = min(r2);

r3 = recall2 - recallPoint;
r3(r3<0) = inf;
[junk,ix2] = min(r3);

PR50 = mean([perc(ix(1)) perc2(ix2(1))]);


%PR50 = PR1;
FP = fp(ix(1));

