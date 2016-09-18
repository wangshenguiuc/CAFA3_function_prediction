

our_pred= load(['..\analysis\CAFA2\github\CAFA2-master\baselines\hpo\our_most_popular'],'pred');
% save(['..\analysis\CAFA2\github\CAFA2-master\baselines\',tp,'\our_most_popular'],'our_pred.pred');
config = cafa_parse_config('example.job');
cafa_driver_preeval(config)
cafa_driver_eval(config)

tmp_eval=load('E:\swang141\project\SequencingNetwork\Sheng\analysis\CAFA2\github\CAFA2-master\seq-centric\hpo\our_most_popular.mat');


our_oa=load(['..\analysis\CAFA2\github\CAFA2-master\benchmark\groundtruth\hpoa.mat']);

[Lia,Locb] = ismember(our_pred.pred.object,our_oa.oa.object);

cafa_eval=load('E:\swang141\project\SequencingNetwork\Sheng\analysis\CAFA2\github\CAFA2-master\evaluation\hpo_all_type1_mode1\hpo_all_type1_mode1\our_most_popular.mat');
ant = our_oa.oa.annotation(Locb,:);
pred = our_pred.pred.score;
pred(:,101) = [];
ant(:,101) = [];
eia = our_oa.oa.eia;
eia(101) = [];
our_eia = pfp_eia(our_pred.pred.ontology.DAG, our_oa.oa.annotation(Locb,:));
% cafa_eia = our_oa.oa.eia;
our_eval = protein_evaluation(pred,ant,eia);

x=tmp_eval.rm.metric{10}(1:5,2)
% tmp_eval.cm_seq.cm(5,10).TN
% tmp_eval.cm_seq.cm(5,10).FP
% tmp_eval.cm_seq.cm(5,10).FN
% tmp_eval.cm_seq.cm(5,10).TP
our_eval.ru_curve_det.R(1:5,10)
cafa_eval.seq_rmcurve.curve(1:5,1:2)
our_eval.ru_curve(1:5,1:2)