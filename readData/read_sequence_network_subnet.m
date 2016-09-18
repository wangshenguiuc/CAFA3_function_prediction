function [Seq_net,Gene_name_new,Gene_name_new_rev] = read_sequence_network_subnet(Gene_name,test_ind,test_sid)
addpath '../data/network/sequence/'
addpath 'util'

GD = [];
for sp=1:length(Gene_name)
    if sp==test_sid
        continue;
    end
    D = Gene_name{sp};
    G = keys(D);
   GD = [GD,G];
end
Gene_name_new = containers.Map(GD, 1:length(GD));
Gene_name_new_rev = containers.Map(1:length(GD), GD);
[e1,e2,wt,~] = textread('sequence_network.txt','%s%s%f%s');
e1=upper(e1);
e2 = upper(e2);
filter = isKey(Gene_name{test_sid},e1)&isKey(Gene_name_new,e2);
e1 = e1(filter);
e2 = e2(filter);
wt = wt(filter);
e1 = cell2mat(values(Gene_name{test_sid},e1));
e2 = cell2mat(values(Gene_name_new,e2));
[~,I]=unique([e1,e2],'rows');
e1 = e1(I);
e2 = e2(I);
wt = wt(I);
ngene1 =  length(Gene_name{test_sid});
ngene2 =  length(GD);
Seq_net = sparse(e1,e2,wt,ngene1,ngene2);
Seq_net = Seq_net(test_ind,:);
end