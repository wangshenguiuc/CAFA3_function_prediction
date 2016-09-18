function [final_score] = blast(Seq_net, Gene_GO_annotation )
%CLUSDCA Summary of this function goes here
%   Detailed explanation goes here

[nnode,nlabel] = size(Gene_GO_annotation);
nnode = size(Seq_net,1);

final_score=zeros(nnode,nlabel);

    go_anot_gene = cell(1,nlabel);
    active_label = [];
    for g=1:nlabel
        go_anot_gene{g} = find(Gene_GO_annotation(:,g)>0); 
        if ~isempty(go_anot_gene{g})
            active_label = [active_label,g];
        end
    end
    fprintf('calculate go assignment finished\n');
    
        for g=active_label
            final_score(:, g) = max(Seq_net(:,go_anot_gene{g}),[],2);
        end


end

