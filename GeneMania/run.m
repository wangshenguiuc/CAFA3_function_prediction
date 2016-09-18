% runGeneMANIA.m   --> this is a sample script for creating two networks
% and ruuning GeneMANIA using these networks and a list of positive genes. 
% 
% %% the input files are
% Zhang_expData.txt which is microarray data 
% sage_expData_sumTagCounts.txt
% %%%%
% %% 1. import the two files: 
% 
% zhang = importdata('Zhang_expData.txt');
% sage = importdata('sage_expData_sumTagCounts.txt');
% in these two data files genes are the rows and features are the columns.

% 2. make two association networks in a cell structure with fields data,
% collabels, rowlabels --> this is done so that we can use the
% combineKernel.m function to combine the networks (to make sure they have
% the same gene order)
% 
% K = 50; % set the number of neigbours 
% 
% network{1}.data = makeAssociationKernel(zhang.data', K);
% network{1}.rowlabels = zhang.textdata(2:end,1);
% network{1}.collabels = network{1}.rowlabels;
% 
% network{2}.data = makeAssociationKernel(sage.data',K);
% network{2}.rowlabels = sage.textdata(2:end,1);
% network{2}.collabels = network{2}.rowlabels;
% 
numNetworks = 6;
nnode = 6311;
nlabel = 4240;
% numNetworks = 6;
% nnode = 16662;
% nlabel = 13708;
%% 3. normalize the networks
for ii = 1:numNetworks
    file_name = ['/home/swang141/research/Bio Network Ontology Prediction/GoPrediction/Data/YeastGraph/noisogo.txt'];
    [g1,g2] = textread(file_name, '%d%d');
    n = nnode;
    A = sparse(g1,g2,true,n,n);
    kernels{ii} = normalizeKernel(A);
end
% 
% 
% 
[g1,g2] = textread('/home/swang141/research/Bio Network Ontology Prediction/GoPrediction/Data/YeastGraph/noisonewAnotationAllNode.txt', '%d%d');
labels = sparse(g1,g2,true,nnode,nlabel);
% labels = labels*2-1;


size(labels)
    

% 7. get area under the AUC for prediction, as well as the scores
pp = randperm(nnode); % this is the permutation index for doing the cross-validation;
nFolds = 3; % number of cross-validation folds
[score, areas, pr, weights] = predictMClassesCV(labels, kernels, nFolds,pp);
fprintf('finished\n');
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




