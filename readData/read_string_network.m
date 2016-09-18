function [net_g2i, net_g2i_rev,Gene_net] = read_string_network(specie)
addpath '../data/network/STRING/'
addpath 'util'
network_list = {'coexpression','cooccurrance','database','experimental','fusion','neighborhood'};

all_node_set = [];
for net_type = network_list
    net_type_str = char(net_type);
   [e1,e2,~,~] = textread([specie,'.',net_type_str],'%s%s%f%s');
   e1 = upper(e1);
   e2 = upper(e2);
   node_set = unique([e1',e2']);
%    fprintf('%d nodes in net %s\n',length(node_set),net_type_str);
   all_node_set = unique([all_node_set,node_set]);
end
fprintf('%d nodes in total\n',length(all_node_set));

net_g2i = containers.Map(all_node_set, 1:length(all_node_set));
net_g2i_rev = containers.Map(1:length(all_node_set),all_node_set);
nnet = length(network_list);
Gene_net = cell(nnet,1);

ct = 1;
for net_type = network_list
    net_type_str = char(net_type);
    net = load_network_text([specie,'.',net_type_str],net_g2i);
    Gene_net{ct} = net;
    ct = ct + 1;
end


end