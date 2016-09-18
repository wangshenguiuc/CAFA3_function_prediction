function protein_eval_res = protein_evaluation(class_score,label,eia,global_label,GO_namespace)
% protein_evaluation(class_score,label,global_label,GO_namespace)
g1 = GO_namespace(:,1);
g2 = GO_namespace(:,2);
[nnode,nlabel] = size(label);
go_type = unique(g2);
assert(length(go_type)==3);
ngo_type = length(go_type);
bias = 1;

npos=zeros(1,nlabel);
st=[3,11,31,101,3];
ed=[10,30,100,300,300];
for i=1:nlabel
    npos(i)=nnz(global_label(:,i));
end
tau = [0.00:0.01:1.00];
ntau = length(tau);
X = zeros(nnode,ntau);
Y = zeros(nnode,ntau);
R = zeros(nnode,ntau);
U = zeros(nnode,ntau);
pr_curve = zeros(ntau,2);
ru_curve = zeros(ntau,2);
protein_eval_res.Smin = zeros(1,ngo_type);
protein_eval_res.fmax = zeros(1,ngo_type);
protein_eval_res.mic_auroc = zeros(ngo_type,5);
protein_eval_res.mac_auroc = zeros(ngo_type,5);
for i=go_type'
    func_t=g1(g2==i);
    
    for j=1:5
        cat=func_t(npos(func_t)>=st(j) & npos(func_t)<=ed(j));
        cat = cat';
        cat = 1:nlabel;
        count = 0;
        mac_auc = 0;
        for k=1:nnode
            lb = label(k,cat);
            sc = class_score(k,cat);
            if min(sc)==max(sc)
                sc = zeros(size(sc));
            else
                sc = (sc-min(sc))/(max(sc)-min(sc));
            end
            [X(k,:),Y(k,:),R(k,:),U(k,:)] =  calPRcurve(sc,lb,tau,eia);
            auc = calcROCarea(sc,lb);
            if ~isnan(auc)
                count = count + 1;
                mac_auc = mac_auc + auc;
            end
        end
        
        s=class_score(:,cat);
        s=s(:);
        l=label(:,cat);
        l=l(:);
        mic_auc = calcROCarea(s,l);
        protein_eval_res.mic_auroc(i,j) = mic_auc;
        protein_eval_res.mac_auroc(i,j) = mac_auc/count;
    end
    
    
    protein_eval_res.pr_curve_det.X = X;
    protein_eval_res.pr_curve_det.Y = Y;
    protein_eval_res.ru_curve_det.R = R;
    protein_eval_res.ru_curve_det.U = U;
    for k=1:ntau
        pr_curve(k, 1) = nanmean(X(:,k), 1);
        pr_curve(k, 2) = nanmean(Y(:,k), 1);
        
        ru_curve(k, 1) = nanmean(R(:,k), 1);
        ru_curve(k, 2) = nanmean(U(:,k), 1);
    end
    f1 = (2*pr_curve(:,1).*pr_curve(:,2)./(pr_curve(:,1)+pr_curve(:,2)));
    protein_eval_res.fmax = max(f1);
    S = sqrt(ru_curve(:, 1).^2 + ru_curve(:, 2).^2);
    protein_eval_res.Smin = min(S);
    protein_eval_res.pr_curve = pr_curve;
    protein_eval_res.ru_curve = ru_curve;
end
end
