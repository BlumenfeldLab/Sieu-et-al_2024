% Edited by Jerry. 10/26/2023
function [datamean,datastd] = ScatterPlotBehavior_231026...
    (data,Title,Xlabels,Ylabels,figsize,SaveDir,alpha,ErrType,Yrange)
    % Function to do a scatter plot with error bar of the data 
    % dataTable: a table with labels corresponding to the categories 
    % title: the title of the plot 
    % yaxe: the yaxis label 
    hFig = SetupFigure_220724(figsize);
    colors = [10,139,148;130,10,40;134,179,82;242,110,48;...
            177,224,123;232,214,16]./255;
    Srange = 10000; 
    for j = 1:length(data(1,:))
        datamean(:,j) = mean(data(:,j),'omitnan');
        if ErrType == 1
            datastd(:,j) = std(data(:,j),'omitnan');
        else
            datastd(:,j)=std(data(:,j),'omitnan')/sqrt(length(data(:,j)));
        end
        % Creating a scatter from each conditions
        x = (1/Srange:1/Srange:length(data(:,j))/Srange)+j;
        scatter(x,data(:,j),'SizeData',35,'MarkerFaceColor',...
            'none','MarkerEdgeColor',colors(j,:));
        hold on;
        ee = errorbar(x(round(length(data(:,j))/2)),datamean(:,j),datastd(:,j));
        ee.Marker = '+';
        ee.MarkerSize = 30;
        ee.MarkerFaceColor = 'k';
        ee.MarkerEdgeColor = 'k';
        ee.Color = 'k';
        ee.LineWidth = 2;
        ee.CapSize = 15;
        hold on;
        LineX(:,j) = j+length(data(:,j))/(2*Srange);
    end
    plot(LineX,datamean,'color',[0.7,0.7,0.7]);
    hold off;
    xlim([0.5 4.5]);
    ylim(Yrange);
    set(gca,'XTick',LineX,'XTickLabel',Xlabels);
    title(Title,'FontSize',17,'FontWeight','bold');
    ylabel(Ylabels,'FontSize',15,'FontWeight','bold'); 
    set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold');
    % export and save fig
    FigName = fullfile(SaveDir,Title); % plot name
    saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
    % significance test
    [pval,tbl,stats] = anova1(data);
    if pval < alpha
        c = multcompare(stats,'display','on','ctype','bonferroni');
    else
        c = [];
    end
    % store significance results 
    DataStat.p = pval; 
    DataStat.table = tbl; 
    DataStat.stats = stats; 
    DataStat.MultComp = c; 
    save(fullfile(SaveDir,[Title '-stat.mat']),'DataStat'); 
end
