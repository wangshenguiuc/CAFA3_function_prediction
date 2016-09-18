function [Seq_net,Gene_name,Gene_name_rev] = read_sequence_network(Gene_name)
addpath '../data/network/sequence/'
addpath 'util'

if ~exist('Gene_name','var')
    %create new Gene dict
    [e1,e2,wt,~] = textread('sequence_network.txt','%s%s%f%s');
    e1 = upper(e1);
    e2 = upper(e2);
    all_node_set = unique(union(e1,e2));    
    Gene_name = containers.Map(all_node_set, 1:length(all_node_set));
    ngene =  length(Gene_name);
    Gene_name_rev = containers.Map(1:length(all_node_set),all_node_set);
    e1 = cell2mat(values(Gene_name,e1));
    e2 = cell2mat(values(Gene_name,e2));
    filter = e1>e2;
    e1 = e1(filter);
    e2 = e2(filter);
    wt = wt(filter);
    Seq_net = sparse(e1,e2,wt,ngene,ngene);
    Seq_net = Seq_net + Seq_net';    
else
    [e1,e2,wt,~] = textread('sequence_network.txt','%s%s%f%s');
       e1 = upper(e1);
    e2 = upper(e2);
    ngene =  length(Gene_name);
    filter = isKey(Gene_name,e1)&isKey(Gene_name,e2);
    e1 = e1(filter);
    e2 = e2(filter);
    wt = wt(filter);
    e1 = cell2mat(values(Gene_name,e1));
    e2 = cell2mat(values(Gene_name,e2));
    filter = e1>e2;
    e1 = e1(filter);
    e2 = e2(filter);
    wt = wt(filter);
    [~,I]=unique([e1,e2],'rows');
    e1 = e1(I);
    e2 = e2(I);
    wt = wt(I);
    Seq_net = sparse(e1,e2,wt,ngene,ngene);
    Seq_net = Seq_net + Seq_net';
end
end