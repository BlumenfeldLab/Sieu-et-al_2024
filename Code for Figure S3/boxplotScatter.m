function boxplotScatter(data, Title, yaxe, withlines, sumstats)
% this is a function to plot a box and a scatter plot of the data 
% data is a table with labels corresponding to the categories 
% title is the title of the plot 
% yaxe is the yaxis label 
% withlines is boolean to indicate if to plot with connecting lines 

randx_scatter = false;
% data = dataTable
% color palette to use for plotting 
colors = [10,139,148;130,10,40;134,179,82;242,110,48;177,224,123;232,214,16]./255;

% Creating a boxplot of the data
varNames = data.Properties.VariableNames;   % variable/category names

% convert the table to an array
data = table2array(data);

nvars = size(data, 2); 

% Getting the handles
a = get(get(gca,'children'),'children');

% Setting the boxes to black
%set(a, 'Color', 'k'); set(a,'LineWidth',2); 
set(gca,'fontsize',15)

% Hodling on to add the scatters
hold on

for i = 1 : nvars
    % Creating a scatter from each conditions
    if randx_scatter
        x(:, i) = (i+(rand(size(data(:,i)))-0.5)/6);
    else
        x(:, i) = ones(size(data(:,i)))*i;
        xlim([1-0.1,nvars+0.1]);
    end
    scatter(gca, x(:, i),...
        data(:,i),'r','filled', 'SizeData', 30, 'MarkerFaceColor',...
        colors(i, :), 'MarkerEdgeColor', 'none');
end

if withlines % 10 Jan 2021 xinyuan edit
%     patchline(x', data', 'FaceColor',[0.5 0.5 0.5], 'LineWidth', 1.3, ...
%         'EdgeColor',[0.5 0.5 0.5], 'EdgeAlpha', 0.2, 'FaceAlpha', 0.2); 
%     patchline() calls patch() which actually draws polygons, not only the lines we need
    for sz_num = 1 : size(x,1)
        for period_num = 1 : size(x,2)-1
            plot([x(sz_num, period_num) x(sz_num, period_num+1)], ...
            [data(sz_num, period_num) data(sz_num,period_num+1)], 'color',[0.7,0.7,0.7]);
        end
    end
end

% make the box plot
if strcmp(sumstats,'median')
    boxplot(data,'Labels',{varNames},'symbol','');
elseif strcmp(sumstats, 'average')
    errorbar([1:length(varNames)],mean(data,'omitnan'),std(data,'omitnan'),'+','LineStyle','none',...
        'CapSize',18,'linewidth', 1,'MarkerSize',10,'color', 'k'); 
    xticks([1:length(varNames)]);
    xticklabels(varNames);
end


% Adding tile of plot and axis
title(Title,'FontSize',17, 'FontWeight','bold');
ylabel(yaxe,'FontSize',15, 'FontWeight','bold'); 

hold off
end