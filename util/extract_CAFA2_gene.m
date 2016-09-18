load('E:\swang141\project\SequencingNetwork\Sheng\analysis\CAFA2\github\CAFA2-master\baselines\hpo\BB4H.mat')

fout = fopen('..\data\CAFA2_target_gene\hpo_target.txt','w');
for i=1:length(pred.object)
    fprintf(fout,'%s\n',char(pred.object(i)));
end
fclose(fout);