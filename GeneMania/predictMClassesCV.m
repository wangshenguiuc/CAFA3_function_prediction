function [areas, pr, weights,loo_sc,test_ind] = predictMClassesCV(labelMat, kernels, nFolds,pp)
%function [areas, pr,  weights] = predictMClassesCV(labelMat, kernels, nFolds,,pp)
%  kernels     -- a cell array of matrices, each matrix is nxn (i.e. size(kernel{1})
%  = [n n])
%  LABELMAT    -- a label matrix of nxd where d is the number of functions (e.g. GO
%  categories)
%  NFOLDS    -- The number of folds of cross-validation to perform
%  PP        -- Optional vector of length length(find(labels)) containing a
%               random permutation of the indices 1:length(find(labels)). 
%               You must specify a weighting scheme in the previous parameter
%               in order to use this argument.
%
% Returns AREAS and Precision at 10% recall, 
% vectors of length NFOLDSxd containing the ROC area on each step
% of cross-validation. weights is the weight of each network (matrix) in kernel
% $Revision: 1.1 $


[NN,d] = size(labelMat); % number of prediction tasks
offset = 0;
recall = 0.1;
loo_sc = inf(NN,d);
if nargin < 4
    pp = randperm(NN);
end
test_ind = cell(1,3);
for cv = 1:1
    labels = labelMat;
    lastElem = min(NN, offset + floor(NN/nFolds));
    ix = pp(offset+1:lastElem);
    offset = lastElem;
    
    labels(ix,:) = 0;
    test_ind{cv} = ix;
    [k, wts] = findKernelWeightsFastSW(labels, kernels);
    nNonzero = 0;
    
    for i = 1:length(k)
        nNonzero = nNonzero + nnz(kernels{k(i)});
    end

    combined = sparse([],[],[],size(labels,1),size(labels,1),nNonzero);

    for i = 1:length(k),
        combined = combined + wts(i) * (kernels{k(i)});
    end;

    combined = sparse(combined);
    
    for mm = 1:d % perform all prediction tasks
        if mod(mm,3000)==0;
            fprintf('finished %f\n',mm/d);
        end
        ll = labels(:,mm)*2-1;
        ll(ix) = 0;
        p = getScoreVectorCG(ll, combined);
		loo_sc(ix,mm) = p(ix);
        areas(cv,mm) = calcROCarea(p(ix), labelMat(ix,mm));
        pr(cv,mm) = calcPR_v2(p(ix),labelMat(ix,mm),recall);
   
    end
    
    
    weights(cv,:) = zeros(1,length(kernels) +1);
    weights(cv,[ (k(:) )']) = wts; 

end
