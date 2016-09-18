function oa = read_blast_annotation(afile,GO_name,GO_net,mode,specie)

[gname,~,goname] = textread(afile,'%s%s%s');

gname = upper(gname);

if mode == 1
is_human =strfind(gname,upper(specie));
xx= 1 - cellfun(@isempty,is_human);

    gname = gname(xx==1);
    goname = goname(xx==1);
    
elseif  mode == 2
    xx = zeros(length(gname),1);
    for specie_nick_l = {'RAT','CAEEL','Drome','HUMAN','Yeast'}
        sp = char(specie_nick_l);
%         if strcmp(sp,specie)
%             continue;
%         end
is_human =strfind(gname,upper(sp));
xx= xx + 1 - cellfun(@isempty,is_human);
    end
    gname = gname(xx>0);
    goname = goname(xx>0);
end

gene_set = unique(gname);

ngene = length(gene_set);
ngo = length(GO_name);
All_Gene_name = containers.Map(gene_set, 1:ngene);
All_Gene_name_rev = containers.Map(1:ngene,gene_set);

oa.object = keys(All_Gene_name);
oa.ontology.DAG = GO_net;
oa.ontology.term = keys(GO_name);

filter = isKey(GO_name,goname);
goname = goname(filter);
gname = gname(filter);
gid = cell2mat(values(All_Gene_name,gname));
goid = cell2mat(values(GO_name,goname));
oa.annotation = logical(sparse(gid,goid,1,ngene,ngo));

old_GO_Gene_mat = oa.annotation';
GO_Gene_mat = oa.annotation';
it = 1;
while it==1 || ~isequal(old_GO_Gene_mat , GO_Gene_mat)
    old_GO_Gene_mat = GO_Gene_mat;
    GO_Gene_mat = GO_net * GO_Gene_mat + old_GO_Gene_mat;
    %     fprintf('it:%d, edge num:%d\n',it,length(find(GO_Gene_mat(:)>0)));
    GO_Gene_mat = double(GO_Gene_mat>=1);
    it = it+1;
end

oa.annotation = GO_Gene_mat';

end

