specie_l = {'Human','Drosophila','Mouse','elegans','Yeast'};
specie_name_l = {'Human','Drosophila','Mouse','C. elegans','Yeast'};
method_l = {'clusDCA','GeneMANIA','blast_clusDCA','our_method'};
figure
metod_n = length(method_l);
for i=1:1
    specie = specie_l{i};
    legend = {'clusDCA','GeneMANIA','Additive','Our method'};
    number = zeros(8,metod_n);
    ebar = zeros(8,metod_n);
    for j=1:metod_n
        tmp = dlmread(['..\result\PSB\',specie,'\',specie,'_',method_l{j},'.txt']);
        number(:,j) = tmp(:,3);
        ebar(:,j) = tmp(:,2)/sqrt(13708);
    end
    grouplabel = {'3-10','11-30','31-100','101-300'};
    specie_name = specie_name_l{i};
    title = [specie_name_l{i},' MF'];
    name = [specie,'_MF_AUC'];
    xlabel = 'Number of annoated genes';
    ylabel = 'AUROC';
    number_t = number(1:4,:);
    ebar_t = ebar(1:4,:);
    plot_bar(title,legend,number_t,ebar_t,xlabel,ylabel,grouplabel,name,0);
    xlabel = 'Number of annoated genes';
    ylabel = 'AUROC';
    title = [specie_name_l{i},' BP'];
    name = [specie,'_BP_AUC'];
    number_t = number(5:8,:);
    ebar_t = ebar(5:8,:);
    plot_bar(title,legend,number_t,ebar_t,xlabel,ylabel,grouplabel,name,0);
end