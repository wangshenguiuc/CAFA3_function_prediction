function [GO_name, GO_name_rev,GO_net,GO_namespace] = read_HPO_network()
addpath '../data/network/core_network'

[e1,~,e2] = textread('hpo_rel.txt','%s%s%s');
node_set = unique([e1',e2']);
net_g2i = containers.Map(node_set, 1:length(node_set));

e1 = cell2mat(values(net_g2i,e1));
e2 = cell2mat(values(net_g2i,e2));
GO_name_rev = containers.Map(1:length(node_set),node_set);
GO_name = net_g2i;
ngo = size(net_g2i,1);
GO_net = sparse(e1,e2,1,ngo,ngo);
GO_namespace = [(1:ngo)',ones(ngo,1)];
end