% Edited by Jerry. 09/21/2022
function ScatterPlotBehavior_220921(dataTable,Title,Yaxe)
    % Function to do a scatter plot with error bar of the data 
    % dataTable: a table with labels corresponding to the categories 
    % title: the title of the plot 
    % yaxe: the yaxis label 
    colors = [10,139,148;130,10,40;134,179,82;242,110,48;...
            177,224,123;232,214,16]./255;
    varNames = categorical(dataTable.Properties.VariableNames);
    data = table2array(dataTable);
    Srange = 10000; 
    for j = 1:length(data(1,:))
        datamean(:,j) = mean(data(:,j),'omitnan');
        datastd(:,j) = std(data(:,j),'omitnan');
        % Creating a scatter from each conditions
        x(:,j) = (1/Srange:1/Srange:length(data(:,j))/Srange)+j;
    end
    % scatter plot
    for j = 1:length(data(1,:))
        LL = length(x(:,j));
        scatter(x(:,j),data(:,j),'SizeData',35,'MarkerFaceColor',...
            'none','MarkerEdgeColor',colors(j,:));
        hold on;
        ee = errorbar(x(round(LL/2),j),datamean(:,j),datastd(:,j));
        ee.Marker = '+';
        ee.MarkerSize = 30;
        ee.MarkerFaceColor = 'k';
        ee.MarkerEdgeColor = 'k';
        ee.Color = 'k';
        ee.LineWidth = 2;
        ee.CapSize = 15;
        hold on;
    end
    % plot transition line
    for j = 1:length(data(1,:)) 
        LineX(:,j) = j+length(data(:,j))/(2*Srange);
    end
    plot(LineX,datamean,'color',[0.7,0.7,0.7]);
    hold off;
    xlim([0.5 4.5]);
    set(gca,'XTick',LineX,'XTickLabel',varNames);
    title(Title,'FontSize',17,'FontWeight','bold');
    ylabel(Yaxe,'FontSize',15,'FontWeight','bold'); 
    set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold'); 
end
