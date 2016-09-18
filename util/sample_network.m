addpath '../data/network/STRING/'
addpath 'util'
network_list = {'coexpression','cooccurrance','database','experimental','fusion','neighborhood'};
specie = 'human';
all_node_set = [];
for prob = [0.2,0.4,0.6,0.8]
    for net_type = network_list
        net_type_str = char(net_type);
        fout = fopen(['../data/network/STRING/sample_network/',specie,'.',net_type_str,'.sample.',num2str(prob)],'w');
        [e1,e2,r1,t1] = textread([specie,'.',net_type_str],'%s%s%f%s');
        for i=1:length(e1)
            if rand(1,1)<prob
                continue;
            end
            fprintf(fout,'%s\t%s\t%f\t%s\n',e1{i},e2{i},r1(i),t1{i});
        end
        fclose(fout);
    end
end
