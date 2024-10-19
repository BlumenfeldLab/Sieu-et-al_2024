% Analysis seizure data (mat) labelled by Anna

% spectrogram/psd plot's units is [dB]: 
% 1. normalization (./) using baseline spectrogram/psd mean (across time)
% 2. mean of normed spectrogram/psd
% 3. 10log(mean)
% 4. plot

% bandpower is scatter-plot of dB-power:
% 1. P = 10log(bandpower)
% 2. mean(P), std(P)/sem(P)
% 3. scatter plot

% AnimalID: Impairedsz, Sparedsz, Allsz, Control
% HCstim: ispi, contra, all


% Edited: Jerry, 7/29/2024
clear; 
close all;
clc;
%% Initialization file saving directory
AnimalID = 'Sparedsz'; % Impairedsz,Sparedsz,Allsz,Control
HCstim = 'ipsi'; % ipsi(HC),contra(HC),all. For control, use all
SaveName = [AnimalID '-' HCstim];
if strcmp(AnimalID,'Control')
    periodLabels = {'Baseline','Control','Post','Late'}; 
else
    periodLabels = {'Baseline','Ictal','PostIctal','Late'}; 
end
Datetime = datestr(now,'dd-mmm-yyyy_HH_MM_SS');
% Labelled sz .mat file directory
DataDir=['E:\DataAnalysis\Anna\' AnimalID]; 
% Directory for saving
TempDir = ['E:\DataAnalysis\JerryLiu\SavedFile\' SaveName];
SaveDir = fullfile(TempDir,Datetime);
if ~exist(SaveDir,'dir')
    mkdir(SaveDir);
end
%% Parameter Initialization
% (HCstim) decide the ('electrode location' + 'sz-endtime'):
% (ipsi)   --> ('ipsi_all'+'ipsi')
% (contra) --> ('contra'+'contra'), 
% (all)    --> ('all'+'max_ipsi_contra').
params = GetParameter_220919(HCstim,SaveDir);
% parse in the relevant parameters 
sps = params.SampleRate; % 1000
PreIctalOffset = params.PreIctalOffset*sps; % baseline 60s*1000
IctalCutoff = params.IctalCutoff*sps; % 15s*1000 ictal cutoff
PostIctalOffset = params.PostIctalOffset*sps; % postictal 30s*1000
RecoveryPeriod = params.RecoveryPeriod*sps; % recovery 60s*1000
ControlOffset = params.ControlOffset*sps; % control 15s*1000
artifactPeriod = params.StimulationArtifactPeriod*sps;% 10s*1000
szTimeSkip = params.SeizureOnsetTimeSkip*sps; % 2.5s*1000
PreStimOffset = params.EPPreStimOffset; 
PostStimOffset = params.EPPostStimOffset; 
% fft number and overlap for spectrogram/psd based on period
nff = 1000; 
nOverlap = params.SpectrogramOverlap;  % 1
%% Get labelled seizure information from .mat files
FileName = ls(fullfile(DataDir,'*.mat')); % .mat file name in char
bilateral = false;
% output three seizure related information (cell)
[szCell,szDir,szTiming] = GetSeizureInfo_220919...
    (AnimalID,FileName,DataDir,params,bilateral,SaveDir);
% convert szTiming from cell to struct for GetBehaviorInfo
rowHeadings = {'AllTimes','AllTimesInds','Lick','Click','ClickInMs',...
    'HC_LFP','LO_LFP','hasWheel','Wheel','hasControl'};
szTiming = cell2struct(szTiming, rowHeadings, 2);
% Get click-based Behavior info
Behavior = GetBehaviorInfo_220920(szCell,szDir,szTiming,params,SaveDir);
%% get length of .mat files
for fn = 1:length(szTiming)
    MatL(fn,:) = szTiming(fn).HC_LFP.length;
end
MatL_mean = mean(MatL)/(60*sps); % [min]
MatL_std = std(MatL)/(60*sps);
MatL_sem = MatL_std/sqrt(length(MatL));
save(fullfile(SaveDir,'SessionlengthMean.mat'),'MatL_mean');
save(fullfile(SaveDir,'SessionlengthStd.mat'),'MatL_std');
%% get sz duration mean+-sd,mean+-sem
szTT = Behavior.BySeizure.IncludedSzDuration;
ictalDuration = szTT(:,2);
szDur_mean = mean(ictalDuration);
szDur_sd = std(ictalDuration);
szDur_sem = szDur_sd/sqrt(length(ictalDuration));
%% Plot Behavior
% parse parameters from the params input struct
firstLickWindow = params.firstLickWindow; 
trialWindow = params.trialWindow; 
maxLickRateWindow = params.maxLickRateWindow;
binSize = params.binSize;
preStimOffset = params.PreStimOffset;
postStimOffset = params.PostStimOffset;
error = 'sem'; % 'sem' or 'std'
nanrm = true; % true: remove nan
% array of time points indicating the start time [s] of each bin
binTimePts = -preStimOffset:binSize:postStimOffset-binSize;
%% Mean lick rate time courses for hits and misses. Fig4-B
figNum = 0;
for hm = 1:2
    if hm == 1  % Hits
        toPlot = Behavior.BySeizure.MeanTimeCourse.Hit./binSize;
        ttl = 'Mean Lick Rate for the Hits';
    else        % Misses
        toPlot = Behavior.BySeizure.MeanTimeCourse.Miss./binSize;
        ttl = 'Mean Lick Rate for the Misses';
    end
    figsize = [9 6]; % fig paper size
    figNum = figNum+1;
    hFig(figNum) = SetupFigure_220724(figsize);
    % plot
    plot(binTimePts,toPlot,'LineWidth',2); hold on;  
    title(ttl);
    ylabel('Mean Lick Rate [licks/sec]');
    xlabel('Time Relative to Water Presentation [s]');
    xlim([binTimePts(1),binTimePts(end)]);
    % add a line to indicate water presentation onset
    xl = xline(0,'--r','Click Onset','LineWidth',1.5); hold off; 
    xl.LabelVerticalAlignment = 'top';
    xl.LabelHorizontalAlignment = 'left';
    xl.LabelOrientation = 'horizontal'; 
    legend(periodLabels);
    set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold');
    % export and save fig
    FigName = fullfile(SaveDir, ttl); % plot name
    saveas(hFig(figNum),FigName,'fig'); saveas(hFig(figNum),FigName,'png');
end
%% Mean Lick Rate Time Courses for Hits and Misses. Fig4-B
% output hit and miss mean and sem, added on 08/22/23
% (same as above but separate figures for each period and include sem)
figNum = 0;
for pd = 1:4
    % parse in the arrays to plot and get their SEM
    hit(:,pd) = Behavior.BySeizure.MeanTimeCourse.Hit(:,pd)/binSize;
    miss(:,pd) = Behavior.BySeizure.MeanTimeCourse.Miss(:,pd)/binSize;
    if strcmp(error,'sem')
        hitSEM(:,pd) = std(Behavior.BySeizure.TimeCourse.Hit...
            (:,:,pd)'./binSize,'omitnan')./...
            sqrt(size(Behavior.BySeizure.TimeCourse.Hit(:,:,pd),2));
        missSEM(:,pd) = std(Behavior.BySeizure.TimeCourse.Miss...
            (:,:,pd)'./binSize,'omitnan')./...
            sqrt(size(Behavior.BySeizure.TimeCourse.Miss(:,:,pd),2));
    else
        hitSEM(:,pd) = std(Behavior.BySeizure.TimeCourse.Hit...
            (:,:,pd)'./binSize,'omitnan');
        missSEM(:,pd) = std(Behavior.BySeizure.TimeCourse.Miss...
            (:,:,pd)'./binSize,'omitnan');
    end
    figNum = figNum+1;
    figsize = [9 6]; % fig paper size
    hFig(figNum) = SetupFigure_220724(figsize);
    % plot the mean time course for the hits
    plot(binTimePts,hit(:,pd),'r','LineWidth',1.7); hold on;
    % plot the mean time course for the misses
    plot(binTimePts,miss(:,pd),'b','LineWidth',1.7); hold on;
    % plot the SEM for the hits
    PlotSem_220706(binTimePts,hit(:,pd),hitSEM(:,pd),gca,'r'); hold on;
    % plot the SEM for the misses
    PlotSem_220706(binTimePts,miss(:,pd),missSEM(:,pd),gca,'b'); hold on;
    % set title and labels
    title(['Hit and Miss Mean Lick Rate During ',periodLabels{pd}]);
    xlabel('Time Relative to Water Presentation [s]')
    ylabel('Mean Lick Rate [licks/s]');
    xlim([binTimePts(1), binTimePts(end)])
    ylim([min(horzcat(hit(:,pd)',miss(:,pd)'))-3, ...
        max(horzcat(hit(:,pd)',miss(:,pd)'))+5]); 
    % show where the click onset is
    xl = xline(0,'--k','Click Onset','LineWidth',1.5); hold off;
    xl.LabelVerticalAlignment = 'top';
    xl.LabelHorizontalAlignment = 'left';
    xl.LabelOrientation = 'horizontal'; 
    legend('Hits','Misses');
    set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold');
    % export and save fig
    FigName = fullfile(SaveDir,...
        ['Hit and Miss Mean Lick Rate During ',periodLabels{pd}]); 
    saveas(hFig(figNum),FigName,'fig'); saveas(hFig(figNum),FigName,'png');
end
%% output Fig 4-B mean+-sem, added on 08/22/2023
% 1, pre-click (-2~-0.1s) data mean and sem, Fig 4-B
% hit pre-click mean and sem
hitpreclk_mean = hit(1:20,:); hitpreclk_sem = hitSEM(1:20,:);
% miss pre-click mean and sem
misspreclk_mean = miss(1:20,:); misspreclk_sem = missSEM(1:20,:);
% save hit mean and sem
save(fullfile(SaveDir,'hitpreclk_mean.mat'),'hitpreclk_mean');
save(fullfile(SaveDir,'hitpreclk_sem.mat'),'hitpreclk_sem');
% save miss mean and sem
save(fullfile(SaveDir,'misspreclk_mean.mat'),'misspreclk_mean');
save(fullfile(SaveDir,'misspreclk_sem.mat'),'misspreclk_sem');
%% for 4 periods mean, hit mean and sem
Periodhitpreclk_mean = mean(hitpreclk_mean,1,'omitnan'); 
Periodhitpreclk_sem = mean(hitpreclk_sem,1,'omitnan');
% for 4 periods mean, miss mean and sem
Periodmisspreclk_mean = mean(misspreclk_mean,1,'omitnan'); 
Periodmisspreclk_sem = mean(misspreclk_sem,1,'omitnan');
% save 4 periods hit mean and sem
save(fullfile(SaveDir,'Periodhitpreclk_mean.mat'),'Periodhitpreclk_mean');
save(fullfile(SaveDir,'Periodhitpreclk_sem.mat'),'Periodhitpreclk_sem');
% save 4 periods miss mean and sem
save(fullfile(SaveDir,'Periodmisspreclk_mean.mat'),'Periodmisspreclk_mean');
save(fullfile(SaveDir,'Periodmisspreclk_sem.mat'),'Periodmisspreclk_sem');
%% unpaired t-test between hit and miss in ictal/postictal/late period
alpha = 0.05;
for i = 1:3 % 'ictal','posticatl','late'
    hitTemp = hitpreclk_mean(:,i+1);
    missTemp = misspreclk_mean(:,i+1);
    [h,p,ci,stats] = ...
        ttest2(hitTemp,missTemp,'Alpha',alpha,'Vartype','unequal');
    % store significance results 
    Ttest2Result.h = h; 
    Ttest2Result.p = p; 
    Ttest2Result.ci = ci;
    Ttest2Result.stats = stats;
    titleTemp = ["ictal","posticatl","late"];
    SaveName = ['Hit-Miss-' char(titleTemp(i)) '-Unpaired T-test'];
    save(fullfile(SaveDir,[SaveName '.mat']),'Ttest2Result');    
end
%% Significance of hit mean lick rate across periods. Fig 4-B
HitMeanLickRate(1,:) = Behavior.BySeizure.MeanTimeCourse.Hit(:,1)/binSize;
HitMeanLickRate(2,:) = Behavior.BySeizure.MeanTimeCourse.Hit(:,2)/binSize;
HitMeanLickRate(3,:) = Behavior.BySeizure.MeanTimeCourse.Hit(:,3)/binSize;
HitMeanLickRate(4,:) = Behavior.BySeizure.MeanTimeCourse.Hit(:,4)/binSize;
alpha = 0.05;
[pval,tbl,stats] = anova1(HitMeanLickRate');
if pval < alpha
    c = multcompare(stats,'display','on','ctype','bonferroni');
else
    c = [];
end
% store significance results 
HitMeanLickR.p = pval; 
HitMeanLickR.table = tbl; 
HitMeanLickR.stats = stats; 
HitMeanLickR.MultComp = c; 
save(fullfile(SaveDir,'HitMeanLickRate-stat.mat'),'HitMeanLickR');
%% Significance of all mean lick rate across periods, after sound, Fig 1-I
alpha = 0.05;
% get the mean lick rate associated with the each period for each seizure
means = squeeze(mean(Behavior.BySeizure.TimeCourse.All(binTimePts >=0 ...
    & binTimePts <= maxLickRateWindow(2),:,1:4),1,'omitnan')./binSize);
% may have to do repeated measures anova?
[pmean,tblmean,statsMean] = anova1(means);
if pmean < alpha
    cmean = multcompare(statsMean,'display','on','ctype','bonferroni');
else
    cmean = [];
end
% store the anova information in the results struct
AllMeanLickRate.p = pmean;
AllMeanLickRate.table = tblmean;
AllMeanLickRate.stats = statsMean;
AllMeanLickRate.MultComp = cmean;
save(fullfile(SaveDir,'AllMeanLickRatePostSound-stat.mat'),'AllMeanLickRate');
%% Pre/Post click lick statistics, Fig 1-I
allTimeCourse = Behavior.BySeizure.TimeCourse.All; 
for pd = 1:4
    % get the time course for this period 
    periodTimeCourse = allTimeCourse(:,:,pd)./binSize; 
    % get the mean and SEM 
    mtemp = mean(periodTimeCourse,2,'omitnan');
    avgPeriodTimeCourse(pd,:) = mtemp;
    if strcmp(error, 'sem')
        stemp = std(periodTimeCourse,0,2,'omitnan')./...
            sqrt(size(periodTimeCourse,2)); 
    else
        stemp = std(periodTimeCourse,0,2,'omitnan');
    end
    semPeriodTimeCourse(pd,:) = stemp;
end
%% 1, pre click (-2~-0.1s) data mean and sem, Fig 1-I
preclk_mean = avgPeriodTimeCourse(:,1:20);
preclk_sem = semPeriodTimeCourse(:,1:20);
% save max lick rate mean and std
save(fullfile(SaveDir,'preclk_mean.mat'),'preclk_mean');
save(fullfile(SaveDir,'preclk_sem.mat'),'preclk_sem');
% added on 08/22/23 for 4 periods
Periodpreclk_mean = mean(preclk_mean,2,'omitnan'); 
Periodpreclk_sem = mean(preclk_sem,2,'omitnan');
%% pre-click (-2~-0.1s) significant test, Fig 1-I
alpha = 0.05;
[pval,tbl,stats] = anova1(preclk_mean');
if pval < alpha
    c = multcompare(stats,'display','off','ctype','bonferroni');
else
    c = [];
end
% store significance results 
PreClkStat.p = pval; 
PreClkStat.table = tbl; 
PreClkStat.stats = stats; 
PreClkStat.MultComp = c; 
save(fullfile(SaveDir,'PreClkStat-stat.mat'),'PreClkStat');
%% 2, post click (0~1s) data mean and sem, Fig 1-I
postclk_mean = avgPeriodTimeCourse(:,21:30);
postclk_sem = semPeriodTimeCourse(:,21:30);
% save max lick rate mean and std
save(fullfile(SaveDir,'postclk_mean.mat'),'postclk_mean');
save(fullfile(SaveDir,'postclk_sem.mat'),'postclk_sem');
%% plot mean lick rate time courses All for each period. Fig 1-I
allTimeCourse = Behavior.BySeizure.TimeCourse.All; 
figNum = 0;
for pd = 1:4
    % get the time course for this period 
    periodTimeCourse = allTimeCourse(:,:,pd)./binSize; 
    % get the mean and SEM 
    avgPeriodTimeCourse = mean(periodTimeCourse,2,'omitnan'); 
    if strcmp(error, 'sem')
        semPeriodTimeCourse = std(periodTimeCourse,0,2,'omitnan')./...
            sqrt(size(periodTimeCourse,2)); 
    else
        semPeriodTimeCourse = std(periodTimeCourse,0,2,'omitnan');
    end
    figNum = figNum+1;
    figsize = [9 6]; % fig paper size
    hFig(figNum) = SetupFigure_220724(figsize);
    % plot 
    plot(binTimePts,avgPeriodTimeCourse,'r','LineWidth',2); hold on;
    % plot the error bars 
    PlotSem_220706(binTimePts,avgPeriodTimeCourse,semPeriodTimeCourse,...
        gca,'b'); hold on;
    % title and labels 
    title(['Mean Lick Rate Timecourse During ',periodLabels{pd}]); 
    xlabel('Time Relative to Water Presentation [s]');
    ylabel('Mean Lick Rate [licks/s]'); 
    xlim([binTimePts(1) binTimePts(end)]); 
    % add a line to indicate click onset 
    xl = xline(0,'--k','Click Onset','LineWidth',1.5); hold off;
    xl.LabelVerticalAlignment = 'top';
    xl.LabelHorizontalAlignment = 'left';
    xl.LabelOrientation = 'horizontal'; 
    set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold');
    % export and save fig
    FigName = fullfile(SaveDir,...
        ['Mean Lick Rate Timecourse During ',periodLabels{pd}]); 
    saveas(hFig(figNum),FigName,'fig'); saveas(hFig(figNum),FigName,'png');
end
%% Scatter plots with error bar during the 4 periods for  Fig 1-J
% (1) max lick rate: this is problemtic since the click interval is
% different in baseline,postical,recovery ('a') and ictal ('d')
Title = 'Maximum Lick Rate'; 
ylabels = 'Max Lick Rate [licks/s]';
figsize = [9 6]; % fig paper size
alpha = 0.05;
ErrType = 0; % 1: mean/std, 0: mean/sem
Yrange = [-2 20];
% max lick rate data
MaxLickRate = Behavior.BySeizure.LickRateAll.Max(:,1:4);
% plot and save
[MLRmean,MLRstd] = ScatterPlotBehavior_231026(MaxLickRate,Title,...
    periodLabels,ylabels,figsize,SaveDir,alpha,ErrType,Yrange);
% save max lick rate mean and std
save(fullfile(SaveDir,'MaxLickRate_mean.mat'),'MLRmean');
if ErrType==1
    save(fullfile(SaveDir,'MaxLickRate_std.mat'),'MLRstd');
else
    save(fullfile(SaveDir,'MaxLickRate_sem.mat'),'MLRstd');
end
%% Scatter plots with error bar during the 4 periods for Fig 1-K
% (2) delay to first lick 
Title = 'Delay to First Lick'; 
ylabels = 'Delay [s]';
figsize = [9 6]; % fig paper size
alpha = 0.05;
ErrType = 0; % 1: mean/std, 0: mean/sem
Yrange = [-1 3];
% (average across all clicks in a period for a seizure)
firstLicksHits = cellfun(@(v)mean(v,'omitnan'),...
    Behavior.ByClick.FirstLickAll.Hit(:, 1:4));
firstLicksMiss = cellfun(@(v)mean(v,'omitnan'),...
    Behavior.ByClick.FirstLickAll.Miss(:, 1:4));
firstLicksHits(firstLicksHits>params.trialWindow) = ...
    params.trialWindow;
firstLicksMiss(firstLicksMiss>params.trialWindow) = ...
    params.trialWindow;
firstLicksAll = cat(3,firstLicksHits,firstLicksMiss);
LickDelay = mean(firstLicksAll,3,'omitnan');
% plot and save
[LDmean,LDstd] = ScatterPlotBehavior_231026(LickDelay,Title,...
    periodLabels,ylabels,figsize,SaveDir,alpha,ErrType,Yrange);
% save max lick rate mean and std
save(fullfile(SaveDir,'LickDelay_mean.mat'),'LDmean');
if ErrType==1
    save(fullfile(SaveDir,'LickDelay_std.mat'),'LDstd');
else
    save(fullfile(SaveDir,'LickDelay_sem.mat'),'LDstd');
end
%% All Wheel speed analysis, Fig 1-L, JerryTest+all
BLthd = 1; % 1mV/ms (10cm/s) baseline run threshold
[AllWheelData,InW,OutW,wheelSopeData] = GetWheelInfo_221114...
    (params,szCell,szDir,szTiming,BLthd,SaveDir);
%% plot mean wheel speed, Fig 1-L
Title = 'Mean Wheel Speed';
Yaxe = 'cm/s';
alpha = 0.05;
figsize = [9 6]; % fig paper size
ErrType = 0; % 1: mean/std, 0: mean/sem
Yrange = [0 40];
[allWheelmean,allWheelstd] = ScatterPlotWheel_231027...
    (AllWheelData(InW,:),Title,Yaxe,figsize,alpha,SaveDir,ErrType,Yrange);
% save all wheel speed mean and std
save(fullfile(SaveDir,'allWheel_mean.mat'),'allWheelmean');
if ErrType==1
    save(fullfile(SaveDir,'allWheel_std.mat'),'allWheelstd');
else
    save(fullfile(SaveDir,'allWheel_sem.mat'),'allWheelstd');
end
%% Control Wheel, use Control data, SupplFig2-F
BLthd = 1; % 1mV/ms (10cm/s) baseline run threshold 
[AllWheelCtrlData,InCtrlW,OutCtrlW,wheelCtrlSopeData] = ...
    GetCtrlWheelInfo_221110...
    (params,szCell,szDir,szTiming,BLthd,SaveDir);
%% plot mean control wheel speed, SupplFig2-F
Title = 'Mean Ctrl Wheel Speed';
Yaxe = 'cm/s';
alpha = 0.05;
figsize = [9 6]; % fig paper size
ErrType = 0; % 1: mean/std, 0: mean/sem
Yrange = [0 30];
[ctrlWheelmean,ctrlWheelstd] = ScatterPlotWheel_231027...
    (AllWheelCtrlData(InCtrlW,:),...
    Title,Yaxe,figsize,alpha,SaveDir,ErrType,Yrange);
% save all wheel speed mean and std
save(fullfile(SaveDir,'ctrlWheel_mean.mat'),'ctrlWheelmean');
if ErrType==1
    save(fullfile(SaveDir,'ctrlWheel_std.mat'),'ctrlWheelstd');
else
    save(fullfile(SaveDir,'ctrlWheel_sem.mat'),'ctrlWheelstd');
end
%% num licks, Edited on 08/10/2023. no all max lick right now.
% 1, mean/max lick number all plot (Fig1)
% 2, mean/max lick number w.r.t hit plot (Fig4)
% 3, mean/max lick number w.r.t miss plot (Fig4)
% get lick numbers responsing to each click in 4 periods 
LicksNumHit = Behavior.ByClick.NumLicksAll.Hit(:,1:4); % cell
LicksNumMiss = Behavior.ByClick.NumLicksAll.Miss(:,1:4); % cell
% lick number mean 
MeanLicksNumHit = cellfun(@(v)mean(v,'omitnan'),LicksNumHit);
MeanLicksNumMiss = cellfun(@(v)mean(v,'omitnan'),LicksNumMiss);
LicksNumAll = cat(3,MeanLicksNumHit,MeanLicksNumMiss); % cell
MeanLicksNumAll = mean(LicksNumAll,3,'omitnan');
% lick number max 
MaxLicksNumHit = cellfun(@(v)max(v),LicksNumHit,'UniformOutput',false);
MaxLicksNumMiss = cellfun(@(v)max(v),LicksNumMiss,'UniformOutput',false);
%% signicance testing for lick number All mean 
alpha = 0.05;
[pval,tbl,stats] = anova1(MeanLicksNumAll);
if pval < alpha
    c = multcompare(stats,'display','off','ctype','bonferroni');
else
    c = [];
end
% store significance results 
AllMeanLickNum.p = pval; 
AllMeanLickNum.table = tbl; 
AllMeanLickNum.stats = stats; 
AllMeanLickNum.MultComp = c; 
save(fullfile(SaveDir,'AllMeanLickNum-stat.mat'),'AllMeanLickNum');
%% Plot AllMeanLickNum 
% make the data into a table
dataTable = array2table(MeanLicksNumAll,'VariableNames',periodLabels);
% remove rows containing NaNs: xinyuan edited Jan 18 
if nanrm 
    dataTable = rmmissing(dataTable);
end
figNum = figNum+1;
figsize = [9 6]; % fig paper size
hFig(figNum) = SetupFigure_220724(figsize);
ScatterPlotBehavior_220921(dataTable,['Mean Number of Licks All'...
    ' Across Periods for All Seizures'],[]);
% export and save fig
FigName = fullfile(SaveDir,'Mean Number of Licks All'); 
saveas(hFig(figNum),FigName,'fig'); saveas(hFig(figNum),FigName,'png');
%% Below is Electrophysiology analysis













%% Below is Electrophysiology analysis
%% Click based spectrogram/psd frequency analysis
% Fig 4-C,D,E
[BClk_SPGM,BClk_PSD,BClk_BP,IClk_SPGM,IClk_PSD,IClk_BP,...
    PClk_SPGM,PClk_PSD,PClk_BP,RClk_SPGM,RClk_PSD,RClk_BP]=...
    AnalyzeClickFrequencyV1_221012...
    (params,Behavior,szCell,szDir,szTiming,SaveDir);
%% Click based spectrogram
% LO, Fig 4-C
BClkLO_SPGM = BClk_SPGM.NBClkLO_SPGM;
HIClkLO_SPGM = IClk_SPGM.NHIClkLO_SPGM;
MIClkLO_SPGM = IClk_SPGM.NMIClkLO_SPGM;
HPClkLO_SPGM = PClk_SPGM.NHPClkLO_SPGM;
MPClkLO_SPGM = PClk_SPGM.NMPClkLO_SPGM;
HRClkLO_SPGM = RClk_SPGM.NHRClkLO_SPGM;
MRClkLO_SPGM = RClk_SPGM.NMRClkLO_SPGM;
B2RClkLO_SPGM = cat(2,BClkLO_SPGM,HIClkLO_SPGM,MIClkLO_SPGM,...
    HPClkLO_SPGM,MPClkLO_SPGM,HRClkLO_SPGM,MRClkLO_SPGM);
% HC
BClkHC_SPGM = BClk_SPGM.NBClkHC_SPGM;
HIClkHC_SPGM = IClk_SPGM.NHIClkHC_SPGM;
MIClkHC_SPGM = IClk_SPGM.NMIClkHC_SPGM;
HPClkHC_SPGM = PClk_SPGM.NHPClkHC_SPGM;
MPClkHC_SPGM = PClk_SPGM.NMPClkHC_SPGM;
HRClkHC_SPGM = RClk_SPGM.NHRClkHC_SPGM;
MRClkHC_SPGM = RClk_SPGM.NMRClkHC_SPGM;
B2RClkHC_SPGM = cat(2,BClkHC_SPGM,HIClkHC_SPGM,MIClkHC_SPGM,...
    HPClkHC_SPGM,MPClkHC_SPGM,HRClkHC_SPGM,MRClkHC_SPGM);
%% Plot LO spectrogram, Fig 4-C
bfy = BClk_SPGM.BClkfy;
btx = BClk_SPGM.BClktx;
freqRange = params.SpectrogramFrequencyRange; % 1~50
Xnames = periodLabels;
figsize = [12 9]; % fig paper size
hFig = SetupFigure_220724(figsize);
% baseline
imagesc(1:4,bfy,10*log10(abs(BClkLO_SPGM)));
colormap(jet);  view(0,90);  hold on;
% ictal
imagesc(6.5:9.5,bfy,10*log10(abs(HIClkLO_SPGM)));
colormap(jet);  view(0,90);  hold on;
imagesc(11:14,bfy,10*log10(abs(MIClkLO_SPGM)));
colormap(jet);  view(0,90);  hold on;
% postictal
imagesc(16.5:19.5,bfy,10*log10(abs(HPClkLO_SPGM)));
colormap(jet);  view(0,90);  hold on;
imagesc(21:24,bfy,10*log10(abs(MPClkLO_SPGM)));
colormap(jet);  view(0,90);  hold on;
% recovery
imagesc(26.5:29.5,bfy,10*log10(abs(HRClkLO_SPGM)));
colormap(jet);  view(0,90);  hold on;
imagesc(31:34,bfy,10*log10(abs(MRClkLO_SPGM)));
colormap(jet);  view(0,90);  hold off;
axis xy;  axis tight; 
ylim(freqRange); % freq range
hcb = colorbar;  caxis([-15 25]);
ylabel(hcb,'Power/Frequency[dB/Hz]','FontSize',14,'FontWeight','bold');
title('Click based LO Sprectrogram'); 
xlabel(' '); ylabel('Frequency[Hz]');
set(gca,'XTick',[2.5 10.25 20.25 30.25],'XTickLabel',Xnames);
set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold');
% export and save fig
FigName = fullfile(SaveDir, 'Click based LO spectrogram'); % plot name
saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
%% Click based psd
% LO, Fig 4-D
BClkLO_PSD = BClk_PSD.NBClkLO_PSD;
HIClkLO_PSD = IClk_PSD.NHIClkLO_PSD;
MIClkLO_PSD = IClk_PSD.NMIClkLO_PSD;
%% plot click based LO psd, Fig 4-D
figsize = [9 6]; % fig paper size
hFig = SetupFigure_220724(figsize);
pName = {'Baseline','Ictal Hit','Ictal Miss'};
freqRange = params.SpectrogramFrequencyRange; % 1~50
powerRange = [-5 15];
bfw = BClk_PSD.BClkfw;
plot(bfw,10*log10(abs(BClkLO_PSD)),'k','LineWidth',2); hold on;
plot(bfw,10*log10(abs(HIClkLO_PSD)),'r','LineWidth',2); hold on;
plot(bfw,10*log10(abs(MIClkLO_PSD)),'g','LineWidth',2); hold off;
xlim(freqRange); % freq range
ylim(powerRange); % power range
title('LO Pwelch PSD Estimate');
xlabel('Frequency [Hz]');
ylabel('Power/Frequency [dB/Hz]');
legend(pName);
set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold');
% export and save fig
FigName = fullfile(SaveDir, 'Click based PSD'); % plot name
saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
%% click based LO bandpower
% LO Delta band, Fig 4-E
ClkLO_DeltaBP(1,:) = BClk_BP.NBClkLO_deltaBP;
ClkLO_DeltaBP(2,:) = IClk_BP.NHIClkLO_deltaBP;
ClkLO_DeltaBP(3,:) = IClk_BP.NMIClkLO_deltaBP;
ClkLO_DeltaBP(4,:) = PClk_BP.NHPClkLO_deltaBP;
ClkLO_DeltaBP(5,:) = PClk_BP.NMPClkLO_deltaBP;
ClkLO_DeltaBP(6,:) = RClk_BP.NHRClkLO_deltaBP;
ClkLO_DeltaBP(7,:) = RClk_BP.NMRClkLO_deltaBP;
% save LO Delta band
save(fullfile(SaveDir,'ClkLO_DeltaBP.mat'),'ClkLO_DeltaBP');
%% mean and std of Late period delta bandpower
for i = 1:length(ClkLO_DeltaBP(:,1))
    y(i,:) = 10*log10(abs(ClkLO_DeltaBP(i,:)));
    RClkLO_deltaBP_mean(i,:) = mean(y(i,:),2,'omitnan');
end
% RHClkLO_deltaBP_mean = mean(ClkLO_DeltaBP(6,:),2,'omitnan');
% RMClkLO_deltaBP_mean = mean(ClkLO_DeltaBP(7,:),2,'omitnan');
%% Unpaired t-test
alpha = 0.05;
[h,p,ci,stats] = ttest2(y(6,:),y(7,:),'Alpha',alpha,'Vartype','unequal');
% store significance results 
Ttest2Result.h = h; 
Ttest2Result.p = p; 
Ttest2Result.ci = ci;
Ttest2Result.stats = stats;
save(fullfile(SaveDir,'LateTtest2Result.mat'),'Ttest2Result');
%% Plot click based LO bandpower
% Fig 4-E, Jerry edited 10/27/2023
% add output option (0/1) for std or sem
% output mean and std/sem
% add ylim
Srange = 10000;
figsize = [9 6]; % fig paper size
Title = 'Click LO-Delta Band Power';
Xnames = periodLabels;
Yname = 'Power/Frequency[dB/Hz]';
ErrType = 0; % 1: mean/std, 0: mean/sem
hFig = SetupFigure_220724(figsize);
for i = 1:length(ClkLO_DeltaBP(:,1))
    y = 10*log10(abs(ClkLO_DeltaBP(i,:)));
    x = (1/Srange:1/Srange:length(y)/Srange)+i;
    clkdatamean(i,:) = mean(y,'omitnan'); % mean or median
    if ErrType == 1
        clkdatastd(i,:) = std(y,'omitnan'); % std
    else
        clkdatastd(i,:) = std(y,'omitnan')/sqrt(length(y)); % sem
    end
    scatter(x,y,'SizeData',35,'MarkerFaceColor','none');
    hold on;
    % errbar
    ee = errorbar(x(:,round(length(y)/2+1)),clkdatamean(i,:),clkdatastd(i,:));
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
xlim([0.5 7.5]);
ylim([-10 20]);
set(gca,'XTick',[1 2.5 4.5 6.5],'XTickLabel',Xnames);
title(Title,'FontSize',17,'FontWeight','bold');
ylabel(Yname,'FontSize',15,'FontWeight','bold'); 
% legend([p1 p2],{'Hits','Misses'});
set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold'); 
% save mean/std
save(fullfile(SaveDir,'clkLODBP_mean.mat'),'clkdatamean');
if ErrType == 1
    save(fullfile(SaveDir,'clkLODBP_std.mat'),'clkdatastd');
else
    save(fullfile(SaveDir,'clkLODBP_sem.mat'),'clkdatastd');
end
% export and save fig
FigName = fullfile(SaveDir, 'Click based LO Delta Bandpower'); % plot name
saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
%% significance test
alpha = 0.01;
Pdata = 10*log10(ClkLO_DeltaBP)';
[pval,tbl,stats] = anova1(Pdata);
if pval < alpha
    c = multcompare(stats,'display','on','ctype','bonferroni');
else
    c = [];
end
% store significance results 
ClkBandPower.p = pval; 
ClkBandPower.table = tbl; 
ClkBandPower.stats = stats; 
ClkBandPower.MultComp = c; 
save(fullfile(SaveDir,[Title(1:4) 'ClkBP-stat.mat']),'ClkBandPower');
%% Periods based spectrogram/psd frequency analysis
% Fig1-E,F,G,H (JerryTest+ipsi/contra)
% Fig3-C,D,E (spared/impared+all)
% SuppFig6 (spared/impared+ipsi/contra)
NormType = 'Unit-Norm'; % 'Unit-Norm' or 'Global-Norm'
if strcmp(NormType,'Unit-Norm')
    [szDur,PeriodLO_SPGM,PeriodLO_PSD,PeriodLO_BP,...
        PeriodHC_SPGM,PeriodHC_PSD,PeriodHC_BP] = ...
        AnalyzePeriodFrequencyV3_220923...
        (params,szCell,szDir,szTiming,SaveDir);
else
    [szDur,PeriodLO_SPGM,PeriodLO_PSD,PeriodLO_BP,...
        PeriodHC_SPGM,PeriodHC_PSD,PeriodHC_BP] = ...
        AnalyzePeriodFrequencyV3_220922...
        (params,szCell,szDir,szTiming,SaveDir);
end
%% scatter plot bandpower
% LO Delta band
LODeltaBP(1,:) = PeriodLO_BP.NormBaselineDeltaBP;
LODeltaBP(2,:) = PeriodLO_BP.NormIctalDeltaBP;
LODeltaBP(3,:) = PeriodLO_BP.NormPostictalDeltaBP;
LODeltaBP(4,:) = PeriodLO_BP.NormRecoveryDeltaBP;
% save LO Delta band
save(fullfile(SaveDir,'LODeltaBP.mat'),'LODeltaBP');
% LO Beta band
LOBetaBP(1,:) = PeriodLO_BP.NormBaselineBetaBP;
LOBetaBP(2,:) = PeriodLO_BP.NormIctalBetaBP;
LOBetaBP(3,:) = PeriodLO_BP.NormPostictalBetaBP;
LOBetaBP(4,:) = PeriodLO_BP.NormRecoveryBetaBP;
% save LO Beta band
save(fullfile(SaveDir,'LOBetaBP.mat'),'LOBetaBP');
% HC Delta band
HCDeltaBP(1,:) = PeriodHC_BP.NormBaselineDeltaBP;
HCDeltaBP(2,:) = PeriodHC_BP.NormIctalDeltaBP;
HCDeltaBP(3,:) = PeriodHC_BP.NormPostictalDeltaBP;
HCDeltaBP(4,:) = PeriodHC_BP.NormRecoveryDeltaBP;
% save HC Delta band
save(fullfile(SaveDir,'HCDeltaBP.mat'),'HCDeltaBP');
% HC Beta band
HCBetaBP(1,:) = PeriodHC_BP.NormBaselineBetaBP;
HCBetaBP(2,:) = PeriodHC_BP.NormIctalBetaBP;
HCBetaBP(3,:) = PeriodHC_BP.NormPostictalBetaBP;
HCBetaBP(4,:) = PeriodHC_BP.NormRecoveryBetaBP;
% save HC Beta band
save(fullfile(SaveDir,'HCBetaBP.mat'),'HCBetaBP');
%% plot Delta band, LO, Fig1-G
alpha = 0.05;
figsize = [9 6]; % fig paper size
Title = 'LO-Delta Band Power';
Yname = 'Power/Frequency[dB/Hz]';
Yrange = [-30 40];
ErrType = 0; % 1: mean/std, 0: mean/sem
[LODBPmean,LODBPstd] = BandPowerPlot_231027...
    (LODeltaBP,Title,Yname,Yrange,figsize,SaveDir,alpha,ErrType);
% save LO Delta band mean and std
save(fullfile(SaveDir,'LODeltaBP_mean.mat'),'LODBPmean');
if ErrType==1
    save(fullfile(SaveDir,'LODeltaBP_std.mat'),'LODBPstd');
else
    save(fullfile(SaveDir,'LODeltaBP_sem.mat'),'LODBPstd');
end
%% plot Delta band, HC, Fig1-G, JerryTest+ipsi, 
% suppl Fig 6, spared/impared + ipsi/contral
alpha = 0.05;
figsize = [9 6]; % fig paper size
Title = 'HC-Delta Band Power';Yname = 'Power/Frequency[dB/Hz]';
Yrange = [-30 40];
ErrType = 0; % 1: mean/std, 0: mean/sem
[HCDBPmean,HCDBPstd] = BandPowerPlot_231027...
    (HCDeltaBP,Title,Yname,Yrange,figsize,SaveDir,alpha,ErrType);
% save LO Delta band mean and std
save(fullfile(SaveDir,'HCDeltaBP_mean.mat'),'HCDBPmean');
if ErrType==1
    save(fullfile(SaveDir,'HCDeltaBP_std.mat'),'HCDBPstd');
else
    save(fullfile(SaveDir,'HCDeltaBP_sem.mat'),'HCDBPstd');
end
%% plot Beta band, LO, Fig1-H, JerryTest+All
figsize = [9 6]; % fig paper size
Title = 'LO-Beta Band Power';
Yname = 'Power/Frequency[dB/Hz]';
Yrange = [-30 40];
ErrType = 0; % 1: mean/std, 0: mean/sem
[LOBBPmean,LOBBPstd] = BandPowerPlot_231027...
    (LOBetaBP,Title,Yname,Yrange,figsize,SaveDir,alpha,ErrType);
% save LO Beta band mean and std
save(fullfile(SaveDir,'LOBetaBP_mean.mat'),'LOBBPmean');
if ErrType==1
    save(fullfile(SaveDir,'LOBetaBP_std.mat'),'LOBBPstd');
else
    save(fullfile(SaveDir,'LOBetaBP_sem.mat'),'LOBBPstd');
end
%% plot Beta band, HC,Fig1-H, suppl Fig 6
Title = 'HC-Beta Band Power';
Yrange = [-30 40];
ErrType = 0; % 1: mean/std, 0: mean/sem
[HCBBPmean,HCBBPstd] = BandPowerPlot_231027...
    (HCBetaBP,Title,Yname,Yrange,figsize,SaveDir,alpha,ErrType);
% save HC Beta band mean and std
save(fullfile(SaveDir,'HCBetaBP_mean.mat'),'HCBBPmean');
if ErrType==1
    save(fullfile(SaveDir,'HCBetaBP_std.mat'),'HCBBPstd');
else
    save(fullfile(SaveDir,'HCBetaBP_sem.mat'),'HCBBPstd');
end
%% szDur plot
edges = 0:1:50;
szDurMean = 0;
figsize = [9 6];  % fig paper size
hFig = SetupFigure_220724(figsize);
if ~isempty(szDur)
    h1 = histogram(vertcat(szDur{:})/1000, edges);
    szDurMean = mean(vertcat(szDur{:})/1000);
else
    h1 = histogram([],edges);
end
ylim([0,30]);
xlabel('Seizure Duration [s]'); ylabel('Counts');
title(['Seizure Duration Distribution ' params.endTimeMethod...
    ' mean ' num2str(szDurMean) 's']);
set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold'); 
% export and save fig
FigName = fullfile(SaveDir, 'Seizure Duration Distribution'); % plot name
saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
%% Get Baseline to recovery spectrogram, Fig1-E
% LO 
BLO_SPGM = PeriodLO_SPGM.NormBaselineLO_SPGM;
ILO_SPGM = PeriodLO_SPGM.NormIctalLO_SPGM;
PLO_SPGM = PeriodLO_SPGM.NormPostictalLO_SPGM;
RLO_SPGM = PeriodLO_SPGM.NormRecoveryLO_SPGM;
B2RLO_SPGM = cat(2,BLO_SPGM,ILO_SPGM,PLO_SPGM,RLO_SPGM);
% HC
BHC_SPGM = PeriodHC_SPGM.NormBaselineHC_SPGM;
IHC_SPGM = PeriodHC_SPGM.NormIctalHC_SPGM;
PHC_SPGM = PeriodHC_SPGM.NormPostictalHC_SPGM;
RHC_SPGM = PeriodHC_SPGM.NormRecoveryHC_SPGM;
B2RHC_SPGM = cat(2,BHC_SPGM,IHC_SPGM,PHC_SPGM,RHC_SPGM);
% LO/HC spectrogram mean
B2RLO_SPGM_mean = mean(B2RLO_SPGM,3,'omitnan');
B2RHC_SPGM_mean = mean(B2RHC_SPGM,3,'omitnan');
%% plot baseline to recovery spectrogram, Fig1-E
freqRange = params.SpectrogramFrequencyRange; % 1~50
bfy = PeriodLO_SPGM.bfy;
b2rtx=(1:(PreIctalOffset+IctalCutoff+PostIctalOffset+RecoveryPeriod)/sps)...
    -PreIctalOffset/sps;
% LO spectrogram
figsize = [12 9]; % fig paper size
hFig = SetupFigure_220724(figsize);
imagesc(b2rtx,bfy,10*log10(abs(B2RLO_SPGM_mean)));
axis xy;  axis tight;  hold on;
xline(0,'-k',{'sz onset'},'FontSize',12,...
    'FontWeight','bold','linewidth',5); hold on;
xline(15,'-k',{'sz cutoff'},'FontSize',12,...
    'FontWeight','bold','linewidth',5); hold on;
xline(45,'-k',{'postical offset'},'FontSize',12,...
    'FontWeight','bold','linewidth',5); hold off;
colormap(jet);  view(0,90);
hcb = colorbar;  caxis([-30 40]);
ylabel(hcb,'Power/Frequency[dB/Hz]','FontSize',14,'FontWeight','bold');
ylim(freqRange); % freq range
title('Baseline to Recovery LO Sprectrogram'); 
xlabel('t[s]');
ylabel('Frequency[Hz]');
set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold'); 
% export and save fig
FigName = fullfile(SaveDir,['PeriodSpectrogram-LO' AnimalID]); 
saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
%% HC spectrogram, Fig1-E
figsize = [12 9]; % fig paper size
hFig = SetupFigure_220724(figsize);
imagesc(b2rtx,bfy,10*log10(abs(B2RHC_SPGM_mean)));
axis xy;  axis tight; 
xline(0,'-k',{'sz onset'},'FontSize',12,...
    'FontWeight','bold','linewidth',5); hold on;
xline(15,'-k',{'sz cutoff'},'FontSize',12,...
    'FontWeight','bold','linewidth',5); hold on;
xline(45,'-k',{'postical offset'},'FontSize',12,...
    'FontWeight','bold','linewidth',5); hold off;
colormap(jet);  view(0,90);
hcb = colorbar;  caxis([-30 40]);
ylabel(hcb,'Power/Frequency[dB/Hz]','FontSize',14,'FontWeight','bold');
ylim(freqRange); % freq range
title('Baseline to Recovery HC Sprectrogram'); 
xlabel('t[s]');
ylabel('Frequency[Hz]');
set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold');
% export and save fig
FigName = fullfile(SaveDir,['PeriodSpectrogram-HC' AnimalID]); 
saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
%% baseline,ictal,postictal,and recovery psd, Fig1-F
% LO
% baseline
BLO_PSD = PeriodLO_PSD.NormBaselineLO_PSD;
BLO_PSD_mean = mean(BLO_PSD,2,'omitnan');
BLO_PSD_median = median(BLO_PSD,2,'omitnan');
BLO_PSD_std = std(BLO_PSD,0,2,'omitnan');
BLO_PSD_pt = prctile(BLO_PSD,[25 75],2);
% ictal
ILO_PSD = PeriodLO_PSD.NormIctalLO_PSD;
ILO_PSD_mean = mean(ILO_PSD,2,'omitnan');
ILO_PSD_median = median(ILO_PSD,2,'omitnan');
ILO_PSD_std = std(ILO_PSD,0,2,'omitnan');
ILO_PSD_pt = prctile(ILO_PSD,[25 75],2);
% postictal
PLO_PSD = PeriodLO_PSD.NormPostictalLO_PSD;
PLO_PSD_mean = mean(PLO_PSD,2,'omitnan');
PLO_PSD_median = median(PLO_PSD,2,'omitnan');
PLO_PSD_std = std(PLO_PSD,0,2,'omitnan');
PLO_PSD_pt = prctile(PLO_PSD,[25 75],2);
% recovery
RLO_PSD = PeriodLO_PSD.NormRecoveryLO_PSD;
RLO_PSD_mean = mean(RLO_PSD,2,'omitnan');
RLO_PSD_median = median(RLO_PSD,2,'omitnan');
RLO_PSD_std = std(RLO_PSD,0,2,'omitnan');
RLO_PSD_pt = prctile(RLO_PSD,[25 75],2);
% HC
% baseline
BHC_PSD = PeriodHC_PSD.NormBaselineHC_PSD;
BHC_PSD_mean = mean(BHC_PSD,2,'omitnan');
BHC_PSD_median = median(BHC_PSD,2,'omitnan');
BHC_PSD_std = std(BHC_PSD,0,2,'omitnan');
BHC_PSD_pt = prctile(BHC_PSD,[25 75],2);
% ictal
IHC_PSD = PeriodHC_PSD.NormIctalHC_PSD;
IHC_PSD_mean = mean(IHC_PSD,2,'omitnan');
IHC_PSD_median = median(IHC_PSD,2,'omitnan');
IHC_PSD_std = std(IHC_PSD,0,2,'omitnan');
IHC_PSD_pt = prctile(IHC_PSD,[25 75],2);
% postictal
PHC_PSD = PeriodHC_PSD.NormPostictalHC_PSD;
PHC_PSD_mean = mean(PHC_PSD,2,'omitnan');
PHC_PSD_median = median(PHC_PSD,2,'omitnan');
PHC_PSD_std = std(PHC_PSD,0,2,'omitnan');
PHC_PSD_pt = prctile(PHC_PSD,[25 75],2);
% recovery
RHC_PSD = PeriodHC_PSD.NormRecoveryHC_PSD;
RHC_PSD_mean = mean(RHC_PSD,2,'omitnan');
RHC_PSD_median = median(RHC_PSD,2,'omitnan');
RHC_PSD_std = std(RHC_PSD,0,2,'omitnan');
RHC_PSD_pt = prctile(RHC_PSD,[25 75],2);
%% plot 4 LO psd mean in one fig, Fig1-F
figsize = [9 6]; % fig paper size
hFig = SetupFigure_220724(figsize);
freqRange = params.SpectrogramFrequencyRange; % 1~50
powerRange = [-5 15];
bfwelch = PeriodLO_PSD.bfwelch;
plot(bfwelch,10*log10(BLO_PSD_mean),'k','LineWidth',2); hold on;
plot(bfwelch,10*log10(ILO_PSD_mean),'r','LineWidth',2); hold on;
plot(bfwelch,10*log10(PLO_PSD_mean),'b','LineWidth',2); hold on;
plot(bfwelch,10*log10(RLO_PSD_mean),'g','LineWidth',2); hold off;
xlim(freqRange); % freq range
ylim(powerRange); % power range
title('LO Pwelch PSD Estimate');
xlabel('Frequency [Hz]');
ylabel('Power/Frequency [dB/Hz]');
legend(periodLabels);
set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold');
% export and save fig
FigName = fullfile(SaveDir, 'LO Pwelch PSD Estimate'); % plot name
saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
%% plot 4 HC psd mean in one fig, Fig1-F
figsize = [9 6]; % fig paper size
hFig = SetupFigure_220724(figsize);
freqRange = params.SpectrogramFrequencyRange; % 1~50
powerRange = [-30 30];
bfwelch = PeriodHC_PSD.bfwelch;
plot(bfwelch,10*log10(BHC_PSD_mean),'k','LineWidth',2); hold on;
plot(bfwelch,10*log10(IHC_PSD_mean),'r','LineWidth',2); hold on;
plot(bfwelch,10*log10(PHC_PSD_mean),'b','LineWidth',2); hold on;
plot(bfwelch,10*log10(RHC_PSD_mean),'g','LineWidth',2); hold off;
xlim(freqRange); % freq range
ylim(powerRange); % power range
title('HC Pwelch PSD Estimate');
xlabel('Frequency [Hz]');
ylabel('Power/Frequency [dB/Hz]');
legend(periodLabels);
set(gca,'FontSize',14,'FontWeight','bold','TitleFontWeight','bold');
% export and save fig
FigName = fullfile(SaveDir, 'HC Pwelch PSD Estimate'); % plot name
saveas(hFig,FigName,'fig'); saveas(hFig,FigName,'png');
%%
close all;
%% Nothing below