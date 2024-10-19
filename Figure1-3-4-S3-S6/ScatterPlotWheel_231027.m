% Fig 1L, scatter plot with error bar of the data
% Edited by Jerry. 10/27/2023
function [datamean,datastd] = ScatterPlotWheel_231027...
    (data,Title,Yaxe,figsize,alpha,SaveDir,ErrType,Yrange)
    % fig setup
    hFig = SetupFigure_220724(figsize);
    colors = [10,139,148;130,10,40;134,179,82;242,110,48;...
            177,224,123;232,214,16]./255;
    Srange = 10000;
    periodLabels = {'Baseline','Ictal','PostIctal','Recovery'};
    data = data*(50/(5000/1000)); % change mV/ms to cm/s, length of wheel 50cm
    for j = 1:length(data(1,:))
        datamean(:,j) = mean(data(:,j),'omitnan');
        if ErrType == 1
            datastd(:,j) = std(data(:,j),'omitnan');
        else
            datastd(:,j)=std(data(:,j),'omitnan')/sqrt(length(data(:,j)));
        end
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
    ylim(Yrange);
    set(gca,'XTick',LineX,'XTickLabel',periodLabels);
    title(Title,'FontSize',17,'FontWeight','bold');
    ylabel(Yaxe,'FontSize',15,'FontWeight','bold'); 
    set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold'); 
    % export and save fig
    FigName = fullfile(SaveDir,Title); % plot name
    saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
    % signicance testing 
    [pval,tbl,stats] = anova1(data);
    if pval < alpha
        c = multcompare(stats,'display','on','ctype','bonferroni');
    else
        c = [];
    end
    % store significance results 
    WheelData.p = pval; 
    WheelData.table = tbl; 
    WheelData.stats = stats; 
    WheelData.MultComp = c; 
    save(fullfile(SaveDir,'WheelDataStats.mat'),'WheelData'); 
end
