function read_result_from_file(result_file)


[model_name,mic1,mic2,mic3,mac1,mac2,mac3] = ...
    textread(result_file,'%s%f%f%f%f%f%f','whitespaces','\t');


name={'MF mic 3-300','BP mic 3-300','CC mic 3-300',...
    'MF mac 3-300','BP mac 3-300','CC mac 3-300'};
data = [mic1,mic2,mic3,mac1,mac2,mac3]';
h=bar(data);
ylim([0.5 0.95])
legend(model_name','location','southoutside');
set(gca,'xticklabel',name)



print('../result/function_prediction/human_sequence.png','-dpng')

end
