function term_eval_res = term_evaluation(class_score,label,global_label,GO_namespace,bias)

g1 = GO_namespace(:,1);
g2 = GO_namespace(:,2);
nlabel = size(label,2);
go_type = unique(g2);
npos=zeros(1,nlabel);
st=[3,11,31,101,3];
ed=[10,30,100,300,300];
for i=1:nlabel
    npos(i)=nnz(global_label(:,i));
end

% for i=1:nnode
%     class_score(i,:) = class_score(i,:) - max(class_score(i,:));
% end
term_eval_res.mic_auroc = zeros(length(go_type),5);
term_eval_res.mac_auroc = zeros(length(go_type),5);
for i=go_type'
    func_t=g1(g2==i);
    for j=1:5
        cat=func_t(npos(func_t)>=st(j) & npos(func_t)<=ed(j));
        cat = cat';
        count = 0;
        %         ap = 0;
        %         apc = 0;
        mac_auc_l=[];
        mac_auc=0;
%         %         ap_l=[];
        for k=cat
            auc = calcROCarea(class_score(:,k),label(:,k));
            if ~isnan(auc)
                count = count + 1;
                mac_auc = mac_auc + auc;
                mac_auc_l=[mac_auc_l,auc];
            end
            if bias==1
                class_score(:,k) = class_score(:,k)-mean(class_score(:,k));
            end

        end
        s=class_score(:,cat);
        s=s(:);
        l=label(:,cat);
        l=l(:);
        mic_auc = calcROCarea(s,l);
        term_eval_res.mic_auroc(i,j) = mic_auc;
        term_eval_res.mac_auroc(i,j) = mac_auc/count;
%         fprintf('type=%d\tcat=%d\tmacro AUC=%f\tstd(macro AUC)=%f\tmicro AUC=%f\tcount=%d\n',i,j,mac_auc/count,std(mac_auc_l),mic_auc,length(cat));
%         
    end
    
end
% 
% mic_auroc = mic_auroc(:,5)';
% mac_auroc = mac_auroc(:,5)';

end
