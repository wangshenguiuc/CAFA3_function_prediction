function [ embedding_vec,net_g2i,net_g2i_rev] = read_embedding_new_Gene(file_name)
%READ_EMBEDDING Summary of this function goes here
%   Detailed explanation goes here
addpath '../data/embedding_vector/LINE/vector'

T = readtable(file_name,'FileType','text','HeaderLines',1,'ReadVariableNames',false,'delimiter',' ');

row_name = table2cell(T(:,1));

net_g2i = containers.Map(row_name, 1:length(row_name));
net_g2i_rev = containers.Map(1:length(row_name),row_name);

ngene = size(row_name,1);

dim = size(T,2) - 1;
if isnan(table2array(T(1,end)))
    dim = dim-1;
end

data = table2array(T(:,2:dim+1));

row_id = cell2mat(values(net_g2i,row_name));


embedding_vec = zeros(ngene,dim);

embedding_vec(row_id,:) = data(:,:);

end

