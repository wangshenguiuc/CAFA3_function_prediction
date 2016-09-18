function [USA] = learn_DCA_vector(network,rspx,dim,net_g2i,net_name)


if iscell(network)
    nnet = size(network,1);
    nnode = size(network{1},1);
    for i=1:nnet
        tA = run_diffusion(network{i}, 'personalized-pagerank', struct('maxiter', 20, 'reset_prob', rspx));
        if i==1
            QA = tA;
            continue
        end
        QA = [QA,tA];
    end
    
    alpha = 1/(nnode);
    QA = log(QA+alpha)-log(alpha);
    
    QA=QA*QA';
    
    fprintf('run SVD d=%d\n',dim);tic
    QA = sparse(QA);
    [U,S] = svds(QA,dim);
    LA = U;
    USA = LA*sqrt(sqrt(S));toc
    
    
    % file_name = [path,'swang141/OutputMatrix/NoIsoStringHumanLX ','d=',num2str(d),' us=',num2str(us),' rsp=',num2str(rspx),'.txt'];
    % dlmwrite(file_name, LA);
        file_name = ['../data/embedding_vector/clusDCA/',char(net_name),'_Gene_','dim',num2str(dim),' rsp_',num2str(rspx),'.U'];
    dlmwrite(file_name, LA); 
    file_name = ['../data/embedding_vector/clusDCA/',char(net_name),'_Gene_','dim',num2str(dim),' rsp_',num2str(rspx),'.US'];
    dlmwrite(file_name, USA); 
    
else
    
    nnode=size(network,1);
    network = network + 0.8*network';
    QA = run_diffusion(network, 'personalized-pagerank', struct('maxiter', 20, 'reset_prob', rspx));
    alpha = 1/(nnode*nnode);
    QA = log(QA+alpha)-log(alpha);
    
    fprintf('run SVD d=%d\n',dim);tic
    [U,S] = svds(QA,dim);
    LA = U;
    USA = LA*sqrt(S);toc
            file_name = ['../data/embedding_vector/clusDCA/GO_','dim',num2str(dim),' rsp_',num2str(rspx),'.U'];
    dlmwrite(file_name, LA); 
        file_name = ['../data/embedding_vector/clusDCA/GO_','dim',num2str(dim),' rsp_',num2str(rspx),'.US'];
    dlmwrite(file_name, USA); 
end

node_names = keys(net_g2i);
node_ids = values(net_g2i,node_names);

net_i2g = containers.Map(node_ids, node_names);

node_id_sorted = values(net_i2g,num2cell(1:nnode))';

T = table(node_id_sorted,USA);
writetable(T,['../data/embedding_vector/clusDCA/',char(net_name),'.embedding'],'Delimiter','\t','WriteVariableNames',false,'FileType','text');

end

