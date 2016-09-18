function [ embedding_vec ] = read_embedding_GO( Gene_name, file_name)
%READ_EMBEDDING Summary of this function goes here
%   Detailed explanation goes here
addpath '../data/embedding_vector/LINE/vector'

T = readtable(file_name,'FileType','text','HeaderLines',0,'ReadVariableNames',false,'delimiter','\t');

row_name = table2cell(T(:,1));

ngene = size(Gene_name,1);

dim = size(T,2) - 1;


filt = Gene_name.isKey(row_name);
data = table2array(T(filt,2:dim+1));

row_name = row_name(filt);

row_id = cell2mat(values(Gene_name,row_name));


embedding_vec = zeros(ngene,dim);

embedding_vec(row_id,:) = data(:,:);

end

