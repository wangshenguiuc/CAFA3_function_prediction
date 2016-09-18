function run_Meng_embedding(Gene_name,Gene_net,specie)
%RUN_MENG_EMBEDDING Summary of this function goes here
%   Detailed explanation goes here
nnet = length(Gene_net);
ngene = size(Gene_net{1},1);
fout = fopen(['../data/embedding_vector/LINE/input/',specie,'.link'],'w');
for i=1:nnet
    cnet = Gene_net{i};
    [gi,gj,s] = find(cnet);
    gid1 = char(values(Gene_name,num2cell(gi)));
    gid2 = char(values(Gene_name,num2cell(gj)));
    nedge = length(s);    
    for k=1:nedge
        if mod(k,100000)==0
            fprintf('finished %f of network %d\n',k/nedge,i);
        end
    fprintf(fout,'%s %s %f %d\n',gid1(k,:),gid2(k,:),s(k),i);
    end
    fprintf('finished writing network:%d\n',i);
end
fclose(fout);

fout = fopen(['../data/embedding_vector/LINE/input/',specie,'.edge'],'w');
for g=1:ngene
    fprintf(fout,'%s %s\n',char(values(Gene_name,num2cell(g))),specie(1));
end
fclose(fout);

end

