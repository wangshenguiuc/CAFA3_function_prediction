function Gene_GO_annotation  = read_annotation(specie,Gene_name,GO_name,GO_net,annot_file)
addpath '../data/annotation/'

if strcmp(specie,'all')
    cur_specie=[{'Drosophila'},{'Elegans'},{'Human'},{'Mouse'},{'Yeast'}];
    GO_Gene_mat = zeros(size(GO_name,1),size(Gene_name,1));
    for s = cur_specie
        [e1,e2,~,~] = textread([char(s),'.GO.selected_code'],'%s%s%f%s');
        e1 = upper(e1);
        e2 = upper(e2);
        filt = GO_name.isKey(e1) & Gene_name.isKey(e2);
        e1 = e1(filt);
        go = cell2mat(values(GO_name,e1));        
        e2 = e2(filt);
        gene = cell2mat(values(Gene_name,e2));
        GO_Gene_mat = GO_Gene_mat + sparse(go,gene,1,size(GO_name,1),size(Gene_name,1));
    end
elseif strcmp(specie,'hpo')
    [e1,e2,~,~] = textread(annot_file,'%s%s%f%s');
    e1 = upper(e1);
    e2 = upper(e2);
    filt = GO_name.isKey(e1) & Gene_name.isKey(e2);
    e1 = e1(filt);
    go = cell2mat(values(GO_name,e1));
    
    e2 = e2(filt);
    gene = cell2mat(values(Gene_name,e2));
    
    GO_Gene_mat = sparse(go,gene,1,size(GO_name,1),size(Gene_name,1));
else
    
    [e1,e2,~,~] = textread([specie,'.GO.selected_code'],'%s%s%f%s');
    
    
    e1 = upper(e1);
    e2 = upper(e2);
    filt = GO_name.isKey(e1) & Gene_name.isKey(e2);
    e1 = e1(filt);
    go = cell2mat(values(GO_name,e1));
    
    e2 = e2(filt);
    gene = cell2mat(values(Gene_name,e2));
    
    GO_Gene_mat = sparse(go,gene,1,size(GO_name,1),size(Gene_name,1));
end

old_GO_Gene_mat = GO_Gene_mat;

it = 1;
while it==1 || ~isequal(old_GO_Gene_mat , GO_Gene_mat)
    old_GO_Gene_mat = GO_Gene_mat;
    GO_Gene_mat = GO_net * GO_Gene_mat + old_GO_Gene_mat;
        fprintf('it:%d, edge num:%d\n',it,length(find(GO_Gene_mat(:)>0)));
    GO_Gene_mat = double(GO_Gene_mat>=1);
    it = it+1;
end

Gene_GO_annotation = GO_Gene_mat';
end