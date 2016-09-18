function plot_bar(title_t,legend_t,number,ebar,xlabel_t,ylabel_t,grouplabel,name,ind)
%PLOT_BAR Summary of this function goes here
%   Detailed explanation goes here
% subplot(4,2,ind);
% axis square;
%cbygmr
rng(2);
H = {'c','y','b','g','r'};
h=bar(number,1,'BarWidth', 1.0);
% error_std = max(rand(size(number))/200,0.0014);
error_std  =ebar;
for i=1:size(number,2)
set(h(i), 'FaceColor', H{i}); 
hold on;
end
% bar(number,'FaceColor',{'c','b','y','g','m','r'});
title(title_t, 'FontSize', 30);
ymin = floor(min(number(:))*10)/10;

ymax = ceil(max(number(:))*20)/20;
ymin = min(ymin,ymax - 0.15);
ylim([ymin,ymax])
set(gca, 'ytick', ymin:0.05:ymax)

% h_legend=legend(legend_t);
% line1=char(legend_t{1},legend_t{2},legend_t{3});
% line2=char(legend_t{4},legend_t{5});
% line={line1,line2};
% 
%    legend_str = []; 
%        for i=1:4, 
%             legend_str = [legend_str; {num2str(i)}]; 
%        end 
%        columnlegend(2, legend_str,'NorthWest'); 
%        
if ind==-1

h_legend=legend(legend_t,'location','north','Orientation','horizontal');
legend boxoff
set(h_legend,'FontSize',20);
end
xlabel(xlabel_t,'fontsize',25);
ylabel(ylabel_t,'fontsize',25);

set(gca, 'XTick', 1:length(grouplabel), 'XTickLabel', grouplabel, 'FontSize', 15);
% set(gca,'yticklabel',strread(num2str(ymin:0.05:1.0),'%s'), 'FontSize', 35);

fig = gcf;
fig.PaperPositionMode = 'auto';
title_t = strrep(title_t, ' ', '_');

numgroups = size(number, 1); 
numbars = size(number, 2); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
      % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
      x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
      errorbar(x, number(:,i), error_std(:,i), 'k', 'linestyle', 'none');
end
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 100 100]); 
print(['../result/PSB/auc/',name],'-dpng','-r0');
% print(['../result/PSB/auc/',name],'-dpdf','-r0');
hold off;
end

