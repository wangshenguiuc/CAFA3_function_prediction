function plot_bar_complementary(title_t,legend_t,number,xlabel_t,ylabel_t,grouplabel,name,ind)
%PLOT_BAR Summary of this function goes here
%   Detailed explanation goes here
% subplot(4,2,ind);
% axis square;
%cbygmr
rng(2);
H = {'g','y','b','c'};
h=bar(number);
for i=1:size(number,2)
set(h(i), 'FaceColor', H{i}); 
hold on;
end
% bar(number,'FaceColor',{'c','b','y','g','m','r'});
title(title_t, 'FontSize', 30);
ymin = floor(min(number(:))*10)/10;
ymax = floor(max(number(:))*10)/10+0.1;
ylim([ymin,ymax])
set(gca, 'ytick', ymin:0.1:ymax)

a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
% Create a vector of '%' signs
     pct = char(ones(size(a,1),1)*'%'); 
% Append the '%' signs after the percentage values
     new_yticks = [char(a),pct];
% 'Reflect the changes on the plot
     set(gca,'yticklabel',new_yticks) 
     


if true
h_legend=legend(legend_t,'location','northoutside','Orientation','horizontal');

set(h_legend,'FontSize',10);
end

% h_legend=legend(legend_t,'location','northwest','Orientation','vertical');
% set(h_legend,'FontSize',14);

legend boxoff 
xlabel(xlabel_t,'fontsize',20);
ylabel(ylabel_t,'fontsize',20);

set(gca, 'XTick', 1:length(grouplabel), 'XTickLabel', grouplabel, 'FontSize', 15);

fig = gcf;
title_t = strrep(title_t, ' ', '_');

% print(['../result/PSB/auc/',name],'-dpng','-r0');
% print(['../result/PSB/auc/',name],'-dpdf','-r0');
hold off;
end

