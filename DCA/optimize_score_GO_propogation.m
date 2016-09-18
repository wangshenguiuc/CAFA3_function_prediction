function [ new_score ] = optimize_score_GO_propogation(score, GO_net)
%OPTIMIZE_SCORE_GO_PROPOGATION Summary of this function goes here
%   Detailed explanation goes here

old_GO_net = GO_net;
it = 1;
while it==1 || ~isequal(old_GO_net , GO_net)
    old_GO_net = GO_net;
    GO_net = GO_net * GO_net + old_GO_net;
    %     fprintf('it:%d, edge num:%d\n',it,length(find(GO_Gene_mat(:)>0)));
    GO_net = double(GO_net>=1);
    it = it+1;
end

fprintf('progation step:%d\n',it);

[nnode,nlabel] = size(score);
new_score = zeros(nnode,nlabel);
GO_net(1:nlabel+1:nlabel*nlabel) = 1;
for i=1:nlabel
    child = GO_net(i,:)~=0;
    new_score(:,i) = max(score(:,child),[],2);
end

end

