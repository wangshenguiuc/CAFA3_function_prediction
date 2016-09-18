function [GO_name, GO_net, net_g2i_rev] = read_CAFA2_GO_network(tp)
addpath '../data/network/core_network/CAFA2_ontology/'

[e1,~,e2] = textread([char(tp),'o_rel.txt'],'%s%s%s');
node_set = unique([e1',e2']);
net_g2i = containers.Map(node_set, 1:length(node_set));
net_g2i_rev = containers.Map(1:length(node_set),node_set);
e1 = cell2mat(values(net_g2i,e1));
e2 = cell2mat(values(net_g2i,e2));

GO_name = net_g2i;
ngo = size(net_g2i,1);
GO_net = sparse(e1,e2,1,ngo,ngo);

end