addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'

specie = 'Yeast';

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
    fprintf('%d nodes in net %s\n',length(node_set),net_type_str);
    all_node_set = unique([all_node_set,node_set]);
end
fprintf('%d nodes in total\n',length(all_node_set));

net_g2i = containers.Map(all_node_set, 1:length(all_node_set));
net_g2i_rev = containers.Map(1:length(all_node_set),all_node_set);
nnet = length(network_list);


net_type = {'experimental'};
net_type_str = char(net_type);
net = load_network_text([specie,'.',net_type_str],net_g2i);
Gene_net = net;

nnode = length(Gene_name);

tA = run_diffusion(Gene_net, 'personalized-pagerank', struct('maxiter', 20, 'reset_prob', 0.5));

alpha = 1/(nnode);
QA = log(tA+alpha)-log(alpha);
fnode = fopen('../data/network/rwr_network/yeast_experiment.node.txt','w');
fedge = fopen('../data/network/rwr_network/yeast_experiment.edge.txt','w');
[x,y,s] = find(QA);
for i=1:length(s)
    fprintf(fedge,'%s\t%s\t%f\tyeast_experiment\n',net_g2i_rev(x(i)),net_g2i_rev(y(i)),s(i));
end
for i=1:nnode
    fprintf(fnode,'%s\tY\n',net_g2i_rev(i));
end
fclose(fedge);
fclose(fnode);
