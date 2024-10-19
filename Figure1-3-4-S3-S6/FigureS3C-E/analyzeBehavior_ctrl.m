function BehaviorResults = analyzeBehavior_ctrl(Behavior, params, varargin)
% this function performs the behavioral analysis on the behavior struct
% outputted by LickTaskBehavior_ctrl and plots the relevant analyses.
% Inputs:
%   - Behavior: the structure outputted by LickTaskBehavior
%   - params: struct with analysis parameters (output of getParams())
%   - figPath: path to save the figures [optional]
% Outputs:
%   - BehaviorResults: a struct containing some relevant analysis results
% Xinyuan Zheng - xinyuan.zheng@yale.edu 
% Abdo Sharaf - abdo.sharaf@yale.edu 
% Date modified: July 21 2022
%**************************************************************************

%% close all figure 
% this is important since we hide the figures at first so we don't want to
% plot on open figures 
close all

%% optional arguments 
global saveFigs; 
saveFigs = false; 
figPath = []; 

if nargin > 2
    saveFigs = true; 
    figPath = varargin{1};
end

%% define and initialize some variables

% parse parameters from the params input struct
firstLickWindow = params.firstLickWindow; 
trialWindow = params.trialWindow; 
maxLickRateWindow = params.maxLickRateWindow;
binSize = params.binSize;
preStimOffset = params.PreStimOffset;
postStimOffset = params.PostStimOffset;

error = 'std';  
error = params.Analysis.BehPlotError;   % 'sem' or 'std'
sumstats = params.Analysis.BehPlotStats; % 'average' or 'median'

narm = true;
withLines = false;
% array of time points indicating the start time (sec) of each bin
binTimePts = -preStimOffset : binSize : postStimOffset - binSize;
% names of periods
periodLabels = {'Baseline', 'Control', 'Post', 'Late'};
% counter for the figures produces
figGroupNum = 1; 

%% Mean hit rate across periods
% % grab a figure handle
% figGroups{figGroupNum} = figure('visible','on');
% 
% % plot average hit rate for each period
% if strcmp(error,'sem')
%     errorbar(nanmean(Behavior.BySeizure.HitRate, 1),...
%         nanstd(Behavior.BySeizure.HitRate, 0, 1)./sqrt(length(Behavior.BySeizure.HitRate)),...
%         'LineWidth', 2);
% else
%     errorbar(nanmean(Behavior.BySeizure.HitRate, 1),...
%         nanstd(Behavior.BySeizure.HitRate, 0, 1),...
%         'LineWidth', 2);
% end
% ylabel('Hit Rate (% Spared Seizures)');
% title('Mean Hit Rate Across Periods');
% 
% % set xtick labels to period names
% set(gca, 'xtick', 1:4, 'xticklabel', periodLabels); xtickangle(45);
% 
% % touch it up
% set(gca, 'FontSize', 14);              % axis font size
% set(gca, 'FontWeight', 'bold');        % axis font weight
% set(gca, 'TitleFontWeight', 'bold');

%% Mean lick rate time courses for hits and misses
figGroupNum = figGroupNum + 1; 
for hm = 1%:2 % hm = 2
    
    if hm == 1  % Hits
        toPlot = Behavior.BySeizure.MeanTimeCourse.Hit./binSize;
        ttl = 'Mean Lick Rate for the Hits';
    else        % Misses
        toPlot = Behavior.BySeizure.MeanTimeCourse.Miss./binSize;
        ttl = 'Mean Lick Rate for the Misses';
    end
    
    % grab a figure handle
    figGroups{figGroupNum}(hm) = figure('visible','on');
    
    % plot
    plot(binTimePts, toPlot,'LineWidth', 2); hold on
    title(ttl);
    ylabel('Mean Lick Rate (licks/sec)');
    xlabel('Time Relative to Water Presentation (sec)');
    xlim([binTimePts(1), binTimePts(end)]);
    
    % add a line to indicate water presentation onset
    xl = xline(0, '--r','Click Onset', 'LineWidth', 1.5);
    xl.LabelVerticalAlignment = 'top';
    xl.LabelHorizontalAlignment = 'left';
    xl.LabelOrientation = 'horizontal'; 
    
    legend(periodLabels);
    
    % touch it up
    set(gca, 'FontSize', 14);              % axis font size
    set(gca, 'FontWeight', 'bold');        % axis font weight
    set(gca, 'TitleFontWeight', 'bold');
    
    hold off 
end

%% Mean Lick Rate Time Courses for Hits and Misses (same as above but separate figures for each period)
% figGroupNum = figGroupNum + 1; 
% for pd = 1:4 % pd = 2
%     % grab a figure handle
%     figGroups{figGroupNum}(pd) = figure('visible','on');
%     
%     % parse in the arrays to plot and get their SEM
%     hit = Behavior.BySeizure.MeanTimeCourse.Hit(:, pd)/binSize;
%     miss = Behavior.BySeizure.MeanTimeCourse.Miss(:, pd)/binSize;
%     
%     if strcmp(error, 'sem')
%         hitSEM = nanstd(Behavior.BySeizure.TimeCourse.Hit(:,:,pd)'./binSize)./...
%             sqrt(size(Behavior.BySeizure.TimeCourse.Hit(:,:,pd), 2));
%         missSEM = nanstd(Behavior.BySeizure.TimeCourse.Miss(:,:,pd)'./binSize)./...
%             sqrt(size(Behavior.BySeizure.TimeCourse.Miss(:,:,pd), 2));
%     else
%         hitSEM = nanstd(Behavior.BySeizure.TimeCourse.Hit(:,:,pd)'./binSize);
%         missSEM = nanstd(Behavior.BySeizure.TimeCourse.Miss(:,:,pd)'./binSize);
%     end
%     
%     SupplFig2.SupplFig2C.(periodLabels{pd}).Mean_lickBeforeStim_Hit = mean(hit(1:20));
%     SupplFig2.SupplFig2C.(periodLabels{pd}).SDAcrossSz_lickBeforeStim_Hit =...
%         std(mean(Behavior.BySeizure.TimeCourse.Hit(1:20,:,pd),1,'omitnan'),'omitnan');
%     SupplFig2.SupplFig2C.(periodLabels{pd}).SEM_lickBeforeStim_Hit =...
%         mean(hitSEM);
%     
%     % plot the mean time course for the hits
%     plot(binTimePts, hit,'r', 'LineWidth', 1.7); hold on
%     % plot the mean time course for the misses
%     plot(binTimePts, miss, 'b', 'LineWidth', 1.7);
%     % plot the SEM for the hits
%     plot_sem(binTimePts, hit, hitSEM, gca, 'r');
%     % plot the SEM for the misses
%     plot_sem(binTimePts, miss, missSEM, gca, 'b');
%     
%     % set title and labels
%     title(['Hit and Miss Mean Lick Rate During ', periodLabels{pd}]);
%     xlabel('Time Relative to Water Presentation (sec)')
%     ylabel('Mean Lick Rate (licks/sec)');
%     xlim([binTimePts(1), binTimePts(end)])
%     ylim([min(horzcat(hit', miss'))-3, max(horzcat(hit', miss'))+5]); 
%     
%     % show where the click onset is
%     xl = xline(0, '--k','Click Onset', 'LineWidth', 1.5);
%     xl.LabelVerticalAlignment = 'top';
%     xl.LabelHorizontalAlignment = 'left';
%     xl.LabelOrientation = 'horizontal'; 
%     
%     % show legend
%     legend('Hits', 'Misses');
%     
%     % touch up
%     set(gca, 'FontSize', 11.5);              % axis font size
%     set(gca, 'FontWeight', 'bold');        % axis font weight
%     set(gca, 'TitleFontWeight', 'bold');
%     
%     hold off
% end

%% Significance of max lick rate across periods
alpha = 0.05;

% get the mean lick rate associated with the each period for each
% seizure
means = squeeze(nanmean(Behavior.BySeizure.TimeCourse.All(binTimePts >=0 ...
    & binTimePts <= maxLickRateWindow(2),:,1:4), 1)./binSize);

% may have to do repeated measures anova?
[pmean, tblmean, statsMean] = anova1(means);
if pmean < alpha
    cmean = multcompare(statsMean, 'display', 'off');
else
    cmean = [];
end

% store the anova information in the results struct
BehaviorResults.Anova.MeanLickRate.p = pmean;
BehaviorResults.Anova.MeanLickRate.table = tblmean;
BehaviorResults.Anova.MeanLickRate.stats = statsMean;
BehaviorResults.Anova.MeanLickRate.MultComp = cmean;

%% Plot duration vs hitrate by seizures and animals 
% figGroupNum = figGroupNum + 1; 
% for i = 1:2
%     % get the data to plot 
%     if i == 1  % by seizure 
%         xplot = Behavior.BySeizure.DurationAndHitRate(:, 3);
%         yplot = Behavior.BySeizure.DurationAndHitRate(:, 4);
%         ttl = 'Effect of Seizure Duration on Hit Rate by Seizure'; 
%     else       % by animal 
%         xplot = Behavior.ByAnimal.DurationAndHitRate(:, 3);
%         yplot = Behavior.ByAnimal.DurationAndHitRate(:, 4);
%         ttl = 'Effect of Seizure Duration on Hit Rate by Animal'; 
%     end
%     
%     % grab a figure handle 
%     figGroups{figGroupNum}(i) = figure('visible','off'); 
%     
%     % plot 
%     plot(xplot, yplot, 'o', 'LineWidth',2); 
%     xlabel('Seizure Duration (sec)'); 
%     ylabel('Hit Rate'); 
%     title(ttl); 
%     xlim([min(xplot), max(xplot)]); 
%     
%     % touch up 
%     set(gca, 'FontSize', 14);              % axis font size
%     set(gca, 'FontWeight', 'bold');        % axis font weight
%     set(gca, 'TitleFontWeight', 'bold');
% end

%% Plot delay to first lick distributions for each period 

% firstLicks = [Behavior.ByClick.FirstLickAll.Hit;...
%     Behavior.ByClick.FirstLickAll.Miss];
% 
% % loop through period and plot distributions 
% figGroupNum = figGroupNum + 1; 
% for pd = 1:4
%     % first licks for this period 
%     periodFirstLicks = firstLicks(:, pd); 
%     periodFirstLicks = vertcat(periodFirstLicks{:}); 
%     % if the first lick happens after the trial window or if no lick
%     % happens at all, replace with the trial window value so we get an idea
%     % of how many such events are there 
%     periodFirstLicks(periodFirstLicks>trialWindow) = trialWindow; 
%     % replace NaNs with the trial window (in case no licks at all)
%     periodFirstLicks(isnan(periodFirstLicks)) = trialWindow; 
%     
%     % grab a figure 
%     figGroups{figGroupNum}(pd) = figure('visible','on'); 
%     
%     % plot the distribution 
%     histogram(periodFirstLicks, 0:0.1:trialWindow); hold on
%     title(['Delay to First Lick During ', periodLabels{pd}]); 
%     xlabel('Delay (sec)'); 
%     ylabel('Count'); 
%     
%     % plot the separator line 
%     xl = xline(firstLickWindow,'--r','Response Window', 'LineWidth', 1.8);
%     xl.LabelVerticalAlignment = 'top';
%     xl.LabelHorizontalAlignment = 'right'; 
%     
%     % touch up 
%     set(gca, 'FontSize', 14);              % axis font size
%     set(gca, 'FontWeight', 'bold');        % axis font weight
%     set(gca, 'TitleFontWeight', 'bold');
%     
%     hold off 
% 
% end

%% plot mean timecourses (across both hits and misses) for each period 
% allTimeCourse = Behavior.BySeizure.TimeCourse.All; 
% 
% % loop through the period and plot timecourses 
% figGroupNum = figGroupNum + 1; 
% for pd = 1:4
%     
%     % get the time course for this period 
%     periodTimeCourse = allTimeCourse(:, :, pd)./binSize; 
%     
%     % get the mean and SEM 
%     avgPeriodTimeCourse = nanmean(periodTimeCourse, 2); 
%     if strcmp(error, 'sem')
%         semPeriodTimeCourse = nanstd(periodTimeCourse, 0, 2)./sqrt(size(periodTimeCourse, 2)); 
%     else
%         semPeriodTimeCourse = nanstd(periodTimeCourse, 0, 2);
%     end
%         
%     % grab a figure handle 
%     figGroups{figGroupNum}(pd) = figure('visible','on'); 
%     
%     % plot 
%     plot(binTimePts, avgPeriodTimeCourse, 'r', 'LineWidth', 2); hold on
%     
%     % plot the error bars 
%     plot_sem(binTimePts, avgPeriodTimeCourse, semPeriodTimeCourse, gca, 'b'); 
%     
%     % title and labels 
%     title(['Mean Lick Rate Timecourse During ', periodLabels{pd}]); 
%     xlabel('Time Relative to Water Presentation (sec)');
%     ylabel('Mean Lick Rate (licks/sec)'); 
%     xlim([binTimePts(1) binTimePts(end)]); 
%     
%     % add a line to indicate click onset 
%     xl = xline(0, '--k','Click Onset', 'LineWidth', 1.5);
%     xl.LabelVerticalAlignment = 'top';
%     xl.LabelHorizontalAlignment = 'left';
%     xl.LabelOrientation = 'horizontal'; 
%     
%     % touch up 
%     set(gca, 'FontSize', 14);              % axis font size
%     set(gca, 'FontWeight', 'bold');        % axis font weight
%     set(gca, 'TitleFontWeight', 'bold');
% 
% end

%% plot mean timecourses (across both hits and misses) for each period 
% loop through the period and plot timecourses 
% figGroupNum = figGroupNum + 1; 
% % grab a figure handle 
% figGroups{figGroupNum}(pd) = figure('visible','on'); 
% 
% for pd = 1:4
%     % get the time course for this period 
%     periodTimeCourse = allTimeCourse(:, :, pd)./binSize; 
%     % get the mean and SEM 
%     avgPeriodTimeCourse = nanmean(periodTimeCourse, 2); 
%     if strcmp(error, 'sem')
%         semPeriodTimeCourse = nanstd(periodTimeCourse, 0, 2)./sqrt(size(periodTimeCourse, 2)); 
%     else
%         semPeriodTimeCourse = nanstd(periodTimeCourse, 0, 2);
%     end
%     h(pd) = ShadedErrorBar(binTimePts, avgPeriodTimeCourse, semPeriodTimeCourse,'-',1); hold on
%     set(h(pd).edge,'LineWidth',1,'LineStyle','-')
%     set(h(pd).mainLine,'LineWidth',2,'LineStyle','-')
%     
%     SupplFig2.SupplFig2C.(periodLabels{pd}).MeanlickBeforeStim_HitMissAll = mean(avgPeriodTimeCourse(1:20));
%     
% end
% 
% % title and labels 
% title(['Mean Lick Rate Timecourse for Control']); 
% xlabel('Time Relative to Water Presentation (sec)');
% ylabel('Mean Lick Rate (licks/sec)'); 
% xlim([binTimePts(1) binTimePts(end)]); 
% 
% % add a line to indicate click onset 
% xl = xline(0, '--k','Click Onset', 'LineWidth', 1.5);
% xl.LabelVerticalAlignment = 'top';
% xl.LabelHorizontalAlignment = 'left';
% xl.LabelOrientation = 'horizontal'; 
% 
% hleg = legend([h(1).mainLine h(2).mainLine h(3).mainLine h(4).mainLine], periodLabels);

%% Plot box scatter plots for the max lick rate, delay to first lick, and num licks during the 4 periods 
ttls = {'Maximum Lick Rate', 'Delay to First Lick', 'Number of Licks', 'Max Estimated Lick Rate (Sliding Window)', ...
    'Average Wheel Speed'}; 
sigFields = {'MaxLickRate', 'Delay', 'NumLicks', 'MaxEstLickRate', 'AvgWheelSpeed'}; 
ylabels = {'Max Lick Rate (licks/sec)', 'Delay (sec)', 'Number of Licks', 'Max Lick Rate (licks/sec)', ...
    'Average Wheel Speed'}; 
for i = 1:3 % i = 1:5
    % i = 2
    if i == 1    % max lick rate
        % get the max lick rate data
        data = Behavior.BySeizure.LickRateAll.Max(:, 1:4);
    elseif i == 2   % delay to first lick
        % get the data for first licks (average across all clicks in a
        % period for a seizure)
        firstLicksHits = cellfun(@nanmean, Behavior.ByClick.FirstLickAll.Hit(:, 1:4));
        firstLicksMiss = cellfun(@nanmean, Behavior.ByClick.FirstLickAll.Miss(:, 1:4));
        firstLicksHits(firstLicksHits > params.trialWindow) = params.trialWindow;
        firstLicksMiss(firstLicksMiss > params.trialWindow) = params.trialWindow;
        firstLicksAll = cat(3, firstLicksHits, firstLicksMiss);
        data = nanmean(firstLicksAll, 3);
    elseif i == 3            % num licks
        % get the data for num licks (average across all clicks in a period
        % for a seizure)
        numLicksHits = cellfun(@nanmean, Behavior.ByClick.NumLicksAll.Hit(:, 1:4));
        numLicksMiss = cellfun(@nanmean, Behavior.ByClick.NumLicksAll.Miss(:, 1:4));
        numLicksAll = cat(3, numLicksHits, numLicksMiss);
        data = nanmean(numLicksAll, 3);
    elseif i == 4   % max estimated lick rate 
        % get the peak of the estimated lick rate timecourse for each click
        maxEstLickRate = cellfun(@max, Behavior.ByClick.TimeCourse.EstimatedLickRate); 
        % initialize the data array 
        data = zeros(size(Behavior.BySeizure.LickRateAll.Max, 1), 4); 
        % loop through files and periods 
        for fl = 1:size(data, 1)
            for pd = 1:4
                data(fl, pd) = nanmean(maxEstLickRate(Behavior.ByClick.ClickTimes(:, 1)==fl & ...
                    Behavior.ByClick.ClickTimes(:, 3) == pd)); 
            end
        end
    else    % average wheel speed 
        whlSpeed = Behavior.ByClick.WheelSpeed.Mean; 
        fls = unique(whlSpeed(:, 1)); 
        data = zeros(length(fls), 4); 
        % loop through files and periods to get averages 
        for fl = 1:size(data, 1)
            for pd = 1:4
                data(fl, pd) = nanmean(whlSpeed(whlSpeed(:, 3)==pd&whlSpeed(:, 1)==fls(fl), 2)); 
            end
        end
    end
    
    % signicance testing 
    [pval, tbl, stats] = anova1(data);
    if pval < alpha
        c = multcompare(stats, 'display', 'off', 'ctype','bonferroni'); % xinyuan edit
    else
        % note: July 2022
        % if anova does not return sig result, there is no point to run multcompare
        % as no groups are significantly different from each other. 
        c = [];
    end
    
    % store significance results 
    BehaviorResults.Anova.(sigFields{i}).p = pval; 
    BehaviorResults.Anova.(sigFields{i}).table = tbl; 
    BehaviorResults.Anova.(sigFields{i}).stats = stats; 
    BehaviorResults.Anova.(sigFields{i}).MultComp = c; 
    
    % make the data into a table
    dataTable = array2table(data,'VariableNames',periodLabels);
    % remove rows containing NaNs: xinyuan edited Jan 18 
    if narm 
        dataTable = rmmissing(dataTable);
    end
    
    % make the box scatter plot: xinyuan edited Jan 12
    figGroupNum = figGroupNum + i; 
    figGroups{figGroupNum} = figure('visible','on');
    boxplotScatter(dataTable, [ttls{i} ' Across Periods for All Seizures'],...
        ylabels{i}, withLines, sumstats);
    
    if i == 1
        for iii_pd = 1:length(periodLabels) % 
            rmmissing_data = table2array(dataTable);
            SupplFig2.SupplFig2D.(periodLabels{iii_pd}).Mean = mean(rmmissing_data(:,iii_pd),'omitnan');
            SupplFig2.SupplFig2D.(periodLabels{iii_pd}).SD = std(rmmissing_data(:,iii_pd),'omitnan');
        end
    elseif i == 2
        for iii_pd = 1:length(periodLabels) % 
            rmmissing_data = table2array(dataTable);
            SupplFig2.SupplFig2E.(periodLabels{iii_pd}).Mean = mean(rmmissing_data(:,iii_pd),'omitnan');
            SupplFig2.SupplFig2E.(periodLabels{iii_pd}).SD = std(rmmissing_data(:,iii_pd),'omitnan');
        end
    end
    
end

%% the distribution of times separating a miss and the next hit 
% arrays to hold the time differences between impaired-impaired and
% impaired-spared, spared-impaired, spared-spared click sequences 
% note that these are just the distributions of click intervals separated
% over categories and periods 
% miss_miss = cell(4, 1); 
% miss_hit = cell(4, 1); 
% hit_hit = cell(4, 1); 
% hit_miss = cell(4, 1); 
% 
% % these will be times between a miss and the next hit 
% miss_next_hit = cell(4, 1); 
% 
% % get the click timing information 
% allClickTimes = Behavior.ByClick.ClickTimes; 
% allClickTimes(allClickTimes(:, 3)==5,:) = [];   % remove the control period 
% 
% % loop through and analyze each period 
% figGroupNum = figGroupNum + 1; 
% for pd = 1:4 
%     
%     pdClicks = allClickTimes(allClickTimes(:, 3)==pd,:); 
%     
%     % get all time differences 
%     allClickDiffs = diff(pdClicks, 1); 
%     
%     % remove between-files differences 
%     allClickDiffs(allClickDiffs(:,1)==1,:) = []; 
%     
%     % get the events for this period
%     % differences for misses followed by hits (this means that the
%     % difference in the 4th column which is 0 for miss, 1 for hit, will be
%     % 1-0=1
%     miss_hit{pd,1} = allClickDiffs(allClickDiffs(:,4)==1,2);   
%     
%     % hits followed by misses
%     hit_miss{pd,1} = allClickDiffs(allClickDiffs(:,4)==-1,2); 
%     
%     % now let's do this for the miss_miss and the hit_hit conditions 
%     allhitIndcs = find(allClickTimes(:,3)==pd & allClickTimes(:,4)==1);  
%     allmissIndcs = find(allClickTimes(:,3)==pd & allClickTimes(:,4)==0);
%     allhitIndcs(allhitIndcs == size(allClickTimes,1)) = [];     % remove if it's the last hit 
%     allmissIndcs(allmissIndcs == size(allClickTimes,1)) = [];   % remove if it's the last miss
%     
%     % get all the events following a hit/miss
%     afterHitTimes = allClickTimes(allhitIndcs + 1,:); 
%     afterMissTimes = allClickTimes(allmissIndcs + 1,:);
%     
%     % compute the time differences 
%     afterHitDiffs = afterHitTimes - allClickTimes(allhitIndcs, :); 
%     afterMissDiffs = afterMissTimes - allClickTimes(allmissIndcs, :); 
%     
%     % retain only the relevant ones 
%     hit_hit{pd, 1} = afterHitDiffs(afterHitTimes(:,end)==1 & afterHitDiffs(:,1)==0,2); % get only the diffs between clicks in the same file 
%     miss_miss{pd, 1} = afterMissDiffs(afterMissTimes(:,end)==0 & afterMissDiffs(:,1)==0,2); 
%     
%     % loop through the click array and get the time intervals between a miss and
%     % the next hit
%     fileNums = unique(pdClicks(:, 1));  % we'll do this for each recording separately 
%     % loop through the files 
%     for file = 1:length(fileNums)
%         % get the clicks for this file and this period 
%         pdfileClicks = pdClicks(pdClicks(:, 1)==fileNums(file), :); 
%         
%         % get the indices of the misses in this period and this file 
%         pdMissIndcs = find(pdfileClicks(:, 4) == 0); % the miss indcs for this period
%         
%         % if no misses skip iteration 
%         if isempty(pdMissIndcs)
%             continue;
%         end
%         
%         % get the reference miss, initially is the first miss 
%         refMiss = pdMissIndcs(1);
%         clk = refMiss;
%         
%         % loop through the clicks and find the time between a miss and the
%         % next hit 
%         while (~isempty(refMiss))
%             % if it's a miss, skip to the next click 
%             if pdfileClicks(clk, 4) == 0 || pdfileClicks(clk, 4) == 2
%                 clk = clk + 1;
%                 if clk > size(pdfileClicks, 1)
%                     break
%                 end
%                 continue;
%                 
%             % if it's a hit, get the interval and update the reference miss
%             elseif pdfileClicks(clk, 4) == 1
%                 timeToNextHit = pdfileClicks(clk, 2) - pdfileClicks(refMiss, 2);
%                 miss_next_hit{pd} = [miss_next_hit{pd}, timeToNextHit];
%                 % update the reference miss, the miss after this hit
%                 refMiss = pdMissIndcs(find(pdMissIndcs > clk, 1));
%                 clk = refMiss;
%             end
%             
%             if clk > size(pdfileClicks, 1)
%                 break;
%             end            
%         end
%     end
%     
%     % grab a figure handle 
%     figGroups{figGroupNum}(pd) = figure('Visible', 'off'); 
%     
%     % plot histogram of miss_next_hit intervals normalized by the average
%     % of click intervals in this period 
%     miss_next_hit_normalized = miss_next_hit{pd} / nanmean(allClickDiffs(:, 2)); 
%     histogram(miss_next_hit_normalized); 
%     
%     % title and labels 
%     title(['Distribution of Time Intervals Separating a Miss and the Next Hit During ',...
%         periodLabels{pd}]);
%     xlabel('Normalized Delay (sec)'); 
%     ylabel('Count'); 
%     
%     % touch up 
%     set(gca, 'FontSize', 14);              % axis font size
%     set(gca, 'FontWeight', 'bold');        % axis font weight
%     set(gca, 'TitleFontWeight', 'bold');
% end
% 
% % plot number of occurences for each click sequence across periods 
% 
% % organize the data which are basically the length of the time difference
% % arrays computed above 
% barData = [cellfun(@length,miss_hit), cellfun(@length,hit_miss),...
%     cellfun(@length, miss_miss), cellfun(@length, hit_hit)]; 
% 
% % now plot 
% figGroupNum = figGroupNum + 1;
% figGroups{figGroupNum} = figure('visible','off'); 
% bar(categorical(periodLabels), barData); 
% legend('Miss-Hit', 'Hit-Miss', 'Miss-Miss', 'Hit-Hit'); 
% title('Number of Occurrences of Consecutive Events'); 
% ylabel('Count'); 
% 
% % touch up 
% set(gca, 'FontSize', 14);              % axis font size
% set(gca, 'FontWeight', 'bold');        % axis font weight
% set(gca, 'TitleFontWeight', 'bold');

%% Plot scatter of delay, lick rate, and number of licks for different periods for all clicks 
% % delay to first lick 
% delay = Behavior.ByClick.FirstLickAll.AllClicks;  
% % total number of licks 
% numLicks = Behavior.ByClick.NumLicksAll.AllClicks;  
% % overall lick rate 
% overallLickRate = Behavior.ByClick.LickRateAll.Overall;   
% % maximum lick rate 
% maxLickRate = Behavior.ByClick.LickRateAll.Max; 
% % maximum wheel speed 
% maxWheelSpeed = Behavior.ByClick.WheelSpeed.Max; 
% % mean wheel speed
meanWheelSpeed = Behavior.ByClick.WheelSpeed.Mean; 
% % estimated lick rate 
% maxEstLickRate = cellfun(@max, Behavior.ByClick.TimeCourse.EstimatedLickRate); 
% % normalized window lick rate 
% normWinRate = Behavior.ByClick.LickRateAll.NormWindowAvg;   
% % remove nans and reset timeout to trial window 
% delay(delay(:, 2)>params.trialWindow, 2) = params.trialWindow; 
% delay(isnan(delay(:, 2)), 2) = params.trialWindow; 
% 
% for an = 1:7
%     
%     if an == 4 && isempty(normWinRate)
%         continue; 
%     end
%     
%     % grab a figure 
%     figGroupNum = figGroupNum + 1;
%     figGroups{figGroupNum} = figure('visible','off');
%     
%     % colors for groups based on classification method 
%     if strcmp(params.ClassificationMethod, 'first_lick')
%         gclrs = ['r', 'g']; 
%         lgnd = {'Miss', 'Hit'}; 
%         mrkSize = 18*ones(1, 2);  
%     else
%         gclrs = ['r', 'g', 'b'];
%         lgnd = {'Miss', 'Hit', 'Unclassified'}; 
%         mrkSize = 18*ones(1, 3); 
%     end
%     
%     % do for each period 
%     for pd = 1:4 
%         subplot(2, 2, pd); 
%         switch an
%             case 1  % plot a scatter of the delay vs overall lick rate 
%                 
%                 bslAvgOverallLickRate = mean(overallLickRate(overallLickRate(:, 3)==pd, 2));
%                 
%                 if strcmp(params.ClassificationMethod, 'norm_lick_rate_and_first_lick') && pd == 1                    
%                     gscatter(delay(delay(:, 3)==pd, 2),...
%                         overallLickRate(overallLickRate(:, 3)==pd, 2)/bslAvgOverallLickRate,...
%                         overallLickRate(overallLickRate(:, 3)==pd, end),...
%                         'g', [], mrkSize);
%                     
%                     legend('Hits');
%                 else
%                     gscatter(delay(delay(:, 3)==pd, 2),...
%                         overallLickRate(overallLickRate(:, 3)==pd, 2)/bslAvgOverallLickRate,...
%                         overallLickRate(overallLickRate(:, 3)==pd, end),...
%                         gclrs, [], mrkSize);
%                     
%                     legend(lgnd);
%                 end
%                                 
%                 % title and labels 
%                 title(['Delay to First Lick vs Normalized Overall Lick Rate During ', periodLabels{pd}]); 
%                 xlabel('Delay (sec)'); 
%                 ylabel('Normalized Overall Lick Rate'); 
%                 
%             case 2  % plot a scatter of the delay vs maximum lick rate 
%                 
%                 if strcmp(params.ClassificationMethod, 'norm_lick_rate_and_first_lick') && pd == 1
%                     gscatter(delay(delay(:, 3)==pd, 2),...
%                         maxLickRate(maxLickRate(:, 3)==pd,2),...
%                         maxLickRate(maxLickRate(:, 3)==pd, end),...
%                         'g', [], mrkSize);
%                     
%                     legend('Hits')
%                     
%                 else
%                     gscatter(delay(delay(:, 3)==pd, 2),...
%                         maxLickRate(maxLickRate(:, 3)==pd,2),...
%                         maxLickRate(maxLickRate(:, 3)==pd, end),...
%                         gclrs, [], mrkSize);
%                     
%                     legend(lgnd)
%                     
%                 end
%                 
%                 % title and labels
%                 title(['Delay to First Lick vs Maximum Lick Rate During ', periodLabels{pd}]); 
%                 xlabel('Delay (sec)'); 
%                 ylabel('Maximum Lick Rate (licks/sec)'); 
%                 
%             case 3  % plot a scatter of the delay vs number of licks 
%                 
%                 if strcmp(params.ClassificationMethod, 'norm_lick_rate_and_first_lick') && pd == 1
%                     gscatter(delay(delay(:, 3)==pd, 2),...
%                         numLicks(numLicks(:, 3)==pd,2),...
%                         numLicks(numLicks(:, 3)==pd, end), 'g', [], mrkSize);
%                     
%                     legend('Hits');
%                     
%                 else
%                     gscatter(delay(delay(:, 3)==pd, 2),...
%                         numLicks(numLicks(:, 3)==pd,2),...
%                         numLicks(numLicks(:, 3)==pd, end), gclrs, [], mrkSize);
%                     
%                     legend(lgnd);
%                     
%                 end
%                 
%                 % title and labels
%                 title(['Delay to First Lick vs Total Number of Licks During ', periodLabels{pd}]);
%                 xlabel('Delay (sec)');
%                 ylabel('Number of Licks');
%                 
%             case 4  % plot a scatter of the delay vs window lick rate 
%                                 
%                 % plot 
%                 if strcmp(params.ClassificationMethod, 'norm_lick_rate_and_first_lick') && pd == 1
%                     gscatter(delay(delay(:, 3)==pd, 2), normWinRate(normWinRate(:, 3)==pd, 2),...
%                         delay(delay(:, 3)==pd, end), 'g', [], mrkSize);
%                     
%                     legend('Hits');
%                     
%                 else
%                     gscatter(delay(delay(:, 3)==pd, 2), normWinRate(normWinRate(:, 3)==pd, 2),...
%                         delay(delay(:, 3)==pd, end), gclrs, [], mrkSize);
%                     
%                     legend(lgnd);
%                     
%                 end
%                 
%                 % title and labels 
%                 title(['Delay to First Lick vs Normalized Window Lick Rate During ', periodLabels{pd}]);
%                 xlabel('Delay (sec)');
%                 ylabel('Window Lick Rate Normalized by Baseline Avg'); 
%                 
%             case 5  % plot delay vs avg estimated lick rate 
%                 
%                 if strcmp(params.ClassificationMethod, 'norm_lick_rate_and_first_lick') && pd == 1
%                     gscatter(delay(delay(:, 3)==pd, 2), maxEstLickRate(delay(:, 3)==pd),...
%                         delay(delay(:, 3)==pd, end), 'g', [], mrkSize);
%                     
%                     legend('Hits')
%                     
%                 else
%                     gscatter(delay(delay(:, 3)==pd, 2), maxEstLickRate(delay(:, 3)==pd),...
%                         delay(delay(:, 3)==pd, end), gclrs, [], mrkSize);
%                     
%                     legend(lgnd)
%                 end
%                 
%                 % title and labels 
%                 title(['Delay to First Lick vs Max Estimated Lick Rate During ', periodLabels{pd}]);
%                 xlabel('Delay (sec)');
%                 ylabel('Max Estimated Lick Rate (licks/sec)'); 
%                 
%             case 6  % plot a scatter of the delay vs max wheel speed 
%                 
%                 if strcmp(params.ClassificationMethod, 'norm_lick_rate_and_first_lick') && pd == 1
%                     gscatter(delay(delay(:, 3)==pd, 2), maxWheelSpeed(delay(:, 3)==pd, 2),...
%                         delay(delay(:, 3)==pd, end), 'g', [], mrkSize);
%                     
%                     legend('Hits')
%                     
%                 else
%                     gscatter(delay(delay(:, 3)==pd, 2), maxWheelSpeed(delay(:, 3)==pd, 2),...
%                         delay(delay(:, 3)==pd, end), gclrs, [], mrkSize);
%                     
%                     legend(lgnd)
%                 end
%                 
%                 % title and labels
%                 title(['Delay to First Lick vs Max Wheel Speed During ', periodLabels{pd}]);
%                 xlabel('Delay (sec)');
%                 ylabel('Max Wheel Speed');
%                 
%             case 7  % plot a scatter of the delay vs mean wheel speed 
%                 
%                 if strcmp(params.ClassificationMethod, 'norm_lick_rate_and_first_lick') && pd == 1
%                     gscatter(delay(delay(:, 3)==pd, 2), meanWheelSpeed(delay(:, 3)==pd, 2),...
%                         delay(delay(:, 3)==pd, end), 'g', [], mrkSize);
%                     
%                     legend('Hits')
%                     
%                 else
%                     gscatter(delay(delay(:, 3)==pd, 2), meanWheelSpeed(delay(:, 3)==pd, 2),...
%                         delay(delay(:, 3)==pd, end), gclrs, [], mrkSize);
%                     
%                     legend(lgnd)
%                 end
%                 
%                 % title and labels
%                 title(['Delay to First Lick vs Mean Wheel Speed During ', periodLabels{pd}]);
%                 xlabel('Delay (sec)');
%                 ylabel('Mean Wheel Speed');
%                  
%         end
%         
%         % touch up
%         set(gca, 'FontSize', 11);              % axis font size
%         set(gca, 'FontWeight', 'bold');        % axis font weight     
%     end               
% end

%% distribution of average wheel speeds for each period 
figGroupNum = figGroupNum + 1;
for pd = 1:4
    % grab a new figure 
    figGroups{figGroupNum}(pd) = figure('visible','off');
    
    % plot histogram 
    histogram(meanWheelSpeed(meanWheelSpeed(:, 3)==pd, 2)); hold on
    title(['Mean Wheel Speed During ', periodLabels{pd}]); 
    xlabel('Mean Wheel Speed'); 
    ylabel('Count'); 
    
    % touch up 
    set(gca, 'FontSize', 14);              % axis font size
    set(gca, 'FontWeight', 'bold');        % axis font weight
    set(gca, 'TitleFontWeight', 'bold');
    
    hold off; 
end

%% standardize figure scales and show 
standardizeScaleAndShow(figGroups, saveFigs, figPath); 
BehaviorResults.SupplFig2 = SupplFig2;

save(fullfile(figPath, 'behResults.mat'),'BehaviorResults');
save(fullfile(figPath, 'SupplFig2.mat'),'SupplFig2');

end

%% Helper functions

function plot_sem(datax, datay, se, ax, clr)
%PLOT_SEM Plot SE/SEM using rectangular drawing function
%   Plot SE/SEM using rectangular drawing function.
%   Currently plots SE; in order to change to SEM just divide SE by the
%   number of points.
%   datax - x values of data on which SEM plot is centered
%   datay - y values of data on which SEM plot is centered
%   se - vector for standard error corresponding to each data
%       point.
%   fig - figure upon which SEM is drawn
%   clr - color (e.g. 'b' - blue, 'g' - green)
if ~isrow(datax)
    datax = datax';
end
if ~isrow(datay)
    datay = datay';
end
if ~isrow(se)
    se = se';
end

X = [datax, flip(datax)];
Y = [datay+se, flip(datay-se)];
fill(ax, X, Y, clr, 'FaceAlpha', 0.1);
end