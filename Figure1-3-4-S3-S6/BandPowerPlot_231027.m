% Band power Plot Fig. 1G,1H

% Edited by Jerry. 10/27/2023
function [datamean,datastd] = BandPowerPlot_231027...
    (BPinput,Title,Yname,Yrange,figsize,SaveDir,alpha,ErrType)
    % fig setup
    hFig = SetupFigure_220724(figsize);
    colors = [10,139,148;130,10,40;134,179,82;242,110,48;...
        177,224,123;232,214,16]./255;
    periodLabels = {'Baseline','Ictal','PostIctal','Recovery'};
    Srange = 10000;
    for i=1:4
        y = 10*log10(BPinput(i,:));
        x = (1/Srange:1/Srange:length(y)/Srange)+i;
        datamean(i,:) = mean(y,'omitnan'); % mean or median
        if ErrType == 1
            datastd(i,:) = std(y,'omitnan'); % std
        else
            datastd(i,:) = std(y,'omitnan')/sqrt(length(y)); % sem
        end
        scatter(x,y,'SizeData',35,'MarkerFaceColor','none',...
            'MarkerEdgeColor',colors(i,:));
        hold on;
        % errbar
        ee = errorbar(x(:,round(length(y)/2+1)),datamean(i,:),datastd(i,:));
        ee.Marker = '+';
        ee.MarkerSize = 15;
        ee.MarkerFaceColor = 'k';
        ee.MarkerEdgeColor = 'k';
        ee.Color = 'k';
        ee.LineWidth = 2;
        ee.CapSize = 15;
        hold on;
    end
    hold off;
    xlim([0.5 4.5]);
    ylim(Yrange);
    for i = 1:4 
        LineX(:,i) = i+length(BPinput(i,:))/(2*Srange);
    end
    set(gca,'XTick',LineX,'XTickLabel',periodLabels);
    title(Title,'FontSize',17,'FontWeight','bold');
    ylabel(Yname,'FontSize',15,'FontWeight','bold'); 
    set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold'); 
    % export and save fig
    FigName = fullfile(SaveDir,Title); % plot name
    saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
    % significance test
    Pdata = 10*log10(BPinput)';
    [pval,tbl,stats] = anova1(Pdata);
    if pval < alpha
        c = multcompare(stats,'display','on','ctype','bonferroni');
    else
        c = [];
    end
    % store significance results 
    BandPower.p = pval; 
    BandPower.table = tbl; 
    BandPower.stats = stats; 
    BandPower.MultComp = c; 
    save(fullfile(SaveDir,[Title(1:4) 'BP-stat.mat']),'BandPower'); 
end