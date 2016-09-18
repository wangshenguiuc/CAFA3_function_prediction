specie_l = {'Yeast','Human','Drosophila','Mouse','elegans','Yeast'};
specie_name_l = {'Yeast','Human','Drosophila','Mouse','C. elegans','Yeast'};
figure
for i=1:6
    specie = specie_l{i};
    mac_auc = dlmread(['..\result\PSB\auc\',specie,'.txt']);
    number = mac_auc(1:4,[7,4,3,5,8]);
    legend = {'Our method','GeneMANIA','clusDCA','BLAST(1080 species)','BLAST(5 species)'};
    grouplabel = {'3-10','11-30','31-10','100-300'};
    specie_name = specie_name_l{i};
    title = [specie_name_l{i},' MF'];
    name = [specie,'_MF_AUC'];
    xlabel = 'Number of annoated genes';
    ylabel = 'AUROC';
    plot_bar(title,legend,number,xlabel,ylabel,grouplabel,name,(i-1)*2+1);
    number = mac_auc(5:8,[7,4,3,5,8]);    
    xlabel = 'Number of annoated genes';
    ylabel = 'AUROC';
    title = [specie_name_l{i},' BP'];
    name = [specie,'_BP_AUC'];
    plot_bar(title,legend,number,xlabel,ylabel,grouplabel,name,(i-1)*2+2);
end