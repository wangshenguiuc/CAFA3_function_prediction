addpath 'DCA/'
addpath 'util/'
addpath 'readData/'
addpath 'LINE/'
addpath 'blast/'
addpath(genpath('SVM/'));
specie = 'Human'; %HMYDE

[GO_name, GO_name_rev,GO_net, GO_namespace] = read_GO_network();
Gene_name = cell(1,5);
Gene_name_rev = cell(1,5);
Gene_net = cell(1,5);
Gene_GO_train_annotation = cell(1,5);
Gene_embedding = cell(1,5);

fname = 'H.ppi_M.ppi_Y.ppi_D.ppi_E.ppi_seq_GOloo_H_M.asc_Y.asc_D.asc_E.asc.node.txt_500_10_0.500000_50_ppi_rwr_1_seq_rwr_1_enum_9_lr_0.025000_mode_1_meta__H7M9G_H7E9G_H7Y9G_H7D9G_H9G8G_H9G';
embedding_file_name =['\\?\E:\swang141\project\SequencingNetwork\Sheng\data\embedding_vector\selected_code\',specie,'\',fname];

for sp=1:5
    spt = char(specie_l(sp));
    [Gene_name{sp},Gene_name_rev{sp}] = read_string_network(spt);
    Gene_GO_train_annotation{sp} = read_annotation(spt,Gene_name{sp},GO_name,GO_net);
    Gene_embedding{sp} = read_embedding(Gene_name{sp}, embedding_file_name);
end


cvrepeat = 3;
gvec = -3:1:0;
cvec = -2:1:2;
logbase = 10;

fprintf('Generating train kernels:\n');
rbfK = cell(length(gvec), 1);
for i = 1:length(gvec)
    fprintf('%d / %d ... ', i, length(gvec)); tic
    rbfK{i} = rbf_kernel(Xtrain, Xtrain, logbase^gvec(i));
    fprintf('done. '); toc
end

% Set up nested CV data
rng(2)

rp = randperm(nnode);
ntest=nnode/nfold;
st = 1;
ed = ntest;
test_ind_nested = rp(st:ed);
train_ind_nested = rp(ed+1:end);


% Nested cross validation
retmax = -inf;
gmax = 1;
cmax = 1;
for gi = 1:length(gvec)
    for ci = 1:length(cvec)
        tt = tic;
        Ktrain = [(1:size(Xtrain, 1))', rbfK{gmax}];
        model = svmtrain(Ytrain, Ktrain, ['-t 4 -q -c ', num2str(logbase^cvec(cmax))]);
        fprintf('Trained SVM on the entire training data. '); toc
        fprintf('Generating test kernel... ');
        Ktest = [(1:size(Xtest, 1))', rbf_kernel(Xtest, Xtrain, logbase^gvec(gmax))];
        [pred, acc, dec] = svmpredict(Ytrain, Ktrain, model, '-q');
        fprintf('done.\n');
    end
end

