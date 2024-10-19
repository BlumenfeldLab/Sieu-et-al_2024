%% Average MUA RMS by period
% xinyuan.zheng@yale.edu

%% paths 
clear; close all;

datapath = '\Data for Figure 2\Behavioral_partial_MUA'; 

svpath_ep = 'F:\xinyuan\Mouse Project\Results';
group = 'MUAGroup';

t = datestr(now,'dd-mmm-yyyy_HH_MM_SS');
if ~exist(fullfile(svpath_ep,[t,'_',group,'_MUA']), 'dir')
   mkdir(fullfile(svpath_ep,[t,'_',group,'_MUA']));
end
svpath_ep = fullfile(svpath_ep,[t,'_',group,'_MUA']);

%% important variables 
mua_fs = 20000; 
lfp_fs = 1000; 
names = string(ls(datapath)); 
names = names(3:end); 
saveFigs = true; 
params = getParams(); 
perAnimal = false; % can be either true or false

%% any desired changes to the params struct 
params.electrodeLocation = 'all'; 
params.endTimeMethod_EP = 'max_ipsi_contra'; 
params.ClassificationMethod = 'first_lick'; 
params.Normalization.Zero = 0; 
params.Analysis.ByAnimalFirst = false; 
params.Analysis.ByAnimalValue = 'average';
params.SampleRate = mua_fs; 
params.OtherSampleRate = lfp_fs; 
errstats = 'sem';

blAdj = 1;

%% get seizure info and seizure paths 
[seizureCell, seizurePaths] = getAllSeizureInfo(names, datapath, params.electrodeLocation, false, svpath_ep);
animals = unique(seizureCell(2:end, 1)); 

%% average MUA 
blTime = params.PreIctalOffset;
analysisTime = 50;
postIctalTime = params.PostIctalOffset;
recoveryTime = params.RecoveryPeriod;

allMUArms1_TOI = nan((blTime+analysisTime)*lfp_fs, length(seizurePaths)); % aligned with szstart
allMUArms1_TOI_end = nan((postIctalTime+recoveryTime)*lfp_fs, length(seizurePaths)); % aligned with szend
allMUADuration = nan(1, length(seizurePaths));

for file = 1:length(seizurePaths)
    
    filepath = seizurePaths{file};
    szInfo = getSeizureTimingInfo_EP(filepath, seizureCell, params);    
    
    szStart = seizureCell{file+1, 3}; 
    szStartIdx_mua = ceil(szStart * mua_fs); % ceil(szInfo.AllTimes(3))
    szEndIdx_mua = ceil(szInfo.AllTimes(4)); 
    
    % Pad MUA signals
    mua_vrms1 = [nan(blTime*mua_fs,1); szInfo.VRMS1.values; nan(analysisTime*mua_fs,1)];
    szStartIdx_mua = szStartIdx_mua + blTime*mua_fs;
    szEndIdx_mua = szEndIdx_mua + blTime*mua_fs;
    
    szDuration_mua = szEndIdx_mua-szStartIdx_mua;
    allMUADuration(file) = szDuration_mua/mua_fs;
    
    % Time of Interest
    mua_vrms1_TOI = nan((blTime+analysisTime)*mua_fs, 1);
    mua_vrms1_TOI(1:blTime*mua_fs+szDuration_mua) = mua_vrms1(szStartIdx_mua - blTime*mua_fs+1: szEndIdx_mua);
    mua_vrms1_TOI_end = mua_vrms1(szEndIdx_mua+1: szEndIdx_mua+(postIctalTime+recoveryTime)*mua_fs);

    % Remove +-1 s around stimulation
    mua_vrms1_TOI((blTime-4)*mua_fs: (blTime+1)*mua_fs,:) = nan;
        
    % Baseline adjustment
    if blAdj
        bl_val = nanmean(mua_vrms1_TOI(1: (blTime-4)*mua_fs));
        mua_vrms1_TOI = mua_vrms1_TOI/bl_val -1;
        mua_vrms1_TOI_end = mua_vrms1_TOI_end/bl_val -1;
    end
    
    % Downsample
    mua_vrms1_TOI_ds = downsample(mua_vrms1_TOI, mua_fs/lfp_fs); 
    mua_vrms1_TOI_end_ds = downsample(mua_vrms1_TOI_end, mua_fs/lfp_fs); 

    % Store in one matrix
    allMUArms1_TOI(:,file) = mua_vrms1_TOI_ds;
    allMUArms1_TOI_end(:,file) = mua_vrms1_TOI_end_ds;

    if 0
        figure; % MUA vRMS 1s
        ShadedErrorBar((1:length(allMUArms1_TOI))./lfp_fs - blTime, nanmean(mua_vrms1_TOI_ds,2),...
                    nanstd(mua_vrms1_TOI_ds,0,2), {'-b', 'LineWidth', 1.3}, 1)
        title('MUA RMS1s Aligned with Seizure Onset');
        ylabel('MUA'); xlabel('Time (s)');
    end
end

%% Scatter Plots
% get average MUA RMS for each period
blMUA_all = allMUArms1_TOI(1:blTime*lfp_fs,:);
blMUA_allmean = mean(blMUA_all,'omitnan');
ictalMUA_all = allMUArms1_TOI(blTime*lfp_fs+1:end,:);
ictalMUA_allmean = mean(ictalMUA_all,'omitnan');
postictalMUA_all = allMUArms1_TOI_end(1:postIctalTime*lfp_fs,:);
postictalMUA_allmean = mean(postictalMUA_all,'omitnan');
recoveryMUA_all = allMUArms1_TOI_end(postIctalTime*lfp_fs+1:end,:);
recoveryMUA_allmean = mean(recoveryMUA_all,'omitnan');

blmean = mean(mean(blMUA_all(1:(blTime-4)*lfp_fs,:),'omitnan'));
ictalmean = mean(mean(ictalMUA_all(1*lfp_fs+1 : 15*lfp_fs,:),'omitnan'));
postmean = mean(mean(postictalMUA_all,'omitnan'));
latemean = mean(mean(recoveryMUA_all,'omitnan'));

if strcmp(errstats,'std')
    blstd = std(mean(blMUA_all(1:(blTime-4)*lfp_fs,:),'omitnan'));
    ictalstd = std(mean(ictalMUA_all(1*lfp_fs+1 : 15*lfp_fs,:),'omitnan'));
    poststd = std(mean(postictalMUA_all,'omitnan'));
    latestd = std(mean(recoveryMUA_all,'omitnan'));

elseif strcmp(errstats,'sem') 
    blstd = std(mean(blMUA_all(1:(blTime-4)*lfp_fs,:),'omitnan'))/sqrt(length(mean(blMUA_all(1:(blTime-4)*lfp_fs,:),'omitnan')));
    ictalstd = std(mean(ictalMUA_all(1*lfp_fs+1 : 15*lfp_fs,:),'omitnan'))/sqrt(length(mean(ictalMUA_all(1*lfp_fs+1 : 15*lfp_fs,:),'omitnan')));
    poststd = std(mean(postictalMUA_all,'omitnan'))/sqrt(length(mean(postictalMUA_all,'omitnan')));
    latestd = std(mean(recoveryMUA_all,'omitnan'))/sqrt(length(mean(recoveryMUA_all,'omitnan')));
end
   
Fig2E_table = array2table([blmean ictalmean postmean latemean; blstd ictalstd poststd latestd],...
    'VariableNames', {'Baseline', 'Ictal','Postictal', 'Recovery'});
writetable(Fig2E_table,['Fig2E_table_',errstats ,'.xlsx'])

% put in a table
periodMUA_allmean(:,1) = blMUA_allmean;
periodMUA_allmean(:,2) = ictalMUA_allmean;
periodMUA_allmean(:,3) = postictalMUA_allmean;
periodMUA_allmean(:,4) = recoveryMUA_allmean;
MUAallmean_table = array2table(periodMUA_allmean, ...
    'VariableNames', {'Baseline', 'Ictal','Postictal', 'Recovery'});

%% plot the scatters
MUA_scatter = figure; 
boxplotScatter(MUAallmean_table, sprintf('Average MUA RMS of Each Animal by Period'),...
   'Baseline-adjusted Mean RMS (percentage change)', true, params.Analysis.MUAstats, errstats);
set(MUA_scatter, 'Position', [0,0,1900,600]);
fig_axs = findall(MUA_scatter,'Type', 'Axes');
savenm = get(fig_axs(1), 'Title'); savenm = regexprep(savenm.String, ' ', '_'); 
saveas(MUA_scatter, fullfile(svpath_ep, [savenm '.png']), 'png');
saveas(MUA_scatter, fullfile(svpath_ep, [savenm '.eps']), 'eps');
saveas(MUA_scatter, fullfile(svpath_ep, [savenm '.fig']), 'fig');

Fig2F_table = array2table([mean(periodMUA_allmean); std(periodMUA_allmean)/sqrt(length(seizurePaths))],...
    'VariableNames', {'Baseline', 'Ictal','Postictal', 'Recovery'});
writetable(Fig2F_table,['Fig2F_table_',errstats,'.xlsx'])

% plot the scatters - no connection lines between periods
MUA_scatter_noline = figure; 
boxplotScatter(MUAallmean_table, sprintf('Average MUA RMS of Each Animal by Period (No lines ver)'),...
   'Baseline-adjusted Mean RMS (percentage change)', false, params.Analysis.MUAstats, errstats);
set(MUA_scatter_noline, 'Position', [0,0,1900,600]);
fig_axs = findall(MUA_scatter_noline,'Type', 'Axes');
savenm = get(fig_axs(1), 'Title'); savenm = regexprep(savenm.String, ' ', '_'); 
saveas(MUA_scatter_noline, fullfile(svpath_ep, [savenm '.png']), 'png');
saveas(MUA_scatter_noline, fullfile(svpath_ep, [savenm '.eps']), 'eps');
saveas(MUA_scatter_noline, fullfile(svpath_ep, [savenm '.fig']), 'fig');

%% Plots
% MUA vRMS 1s
MUA_RMS1 = figure; 
std_allMUArms1_TOI = nanstd(allMUArms1_TOI,0,2);
sem_allMUArms1_TOI = std_allMUArms1_TOI/sqrt(length(seizurePaths));

if strcmp(errstats, 'std')
ShadedErrorBar((0:length(allMUArms1_TOI)-1)./lfp_fs - blTime, nanmean(allMUArms1_TOI,2),...
           std_allMUArms1_TOI,{'-b', 'LineWidth', 1.3}, 1);
elseif strcmp(errstats, 'sem')
    ShadedErrorBar((0:length(allMUArms1_TOI)-1)./lfp_fs - blTime, nanmean(allMUArms1_TOI,2),...
           sem_allMUArms1_TOI,{'-b', 'LineWidth', 1.3}, 1);
end

ylim([-0.5,0.5]);
if blAdj
title(['Baseline Adjusted Average MUA RMS(1s) Aligned with Seizure Onset - Baseline and Ictal (mean duration '...
    num2str(mean(allMUADuration)) ')']); 
ylabel('MUA RMS(1s) Percentage Change');
else
title('Average MUA RMS(1s) Aligned with Seizure Onset');
ylabel('MUA RMS(1s)');
end
xlabel('Time (s)'); xlim([-60,50]); 
xline(0,'-','LineWidth', 2);
xline(-4,'-','-4s','LabelHorizontalAlignment','left','LineWidth', 2);
xline(1,'-','1s','LabelHorizontalAlignment','right','LineWidth', 2);
xline(mean(allMUADuration),'--',['Mean Seizure Duration ' num2str(mean(allMUADuration)) 's'],'LabelHorizontalAlignment','right');
set(MUA_RMS1, 'Position', [0,0,1900,600]);
fig_axs = findall(MUA_RMS1,'Type', 'Axes');
savenm = get(fig_axs(1), 'Title'); savenm = regexprep(savenm.String, ' ', '_'); 
saveas(MUA_RMS1, fullfile(svpath_ep, [savenm '.png']), 'png');
saveas(MUA_RMS1, fullfile(svpath_ep, [savenm '.eps']), 'eps');
saveas(MUA_RMS1, fullfile(svpath_ep, [savenm '.fig']), 'fig');

% aligned with szend
MUA_RMS1_end = figure; 
std_allMUArms1_TOI_end = nanstd(allMUArms1_TOI_end,0,2);
sem_allMUArms1_TOI_end = std_allMUArms1_TOI_end/sqrt(length(seizurePaths));

if strcmp(errstats, 'std')
ShadedErrorBar((0:length(allMUArms1_TOI_end)-1)./lfp_fs, nanmean(allMUArms1_TOI_end,2),...
            nanstd(allMUArms1_TOI_end,0,2), {'-b', 'LineWidth', 1.3}, 1);
elseif strcmp(errstats, 'sem')
ShadedErrorBar((0:length(allMUArms1_TOI_end)-1)./lfp_fs, nanmean(allMUArms1_TOI_end,2),...
            sem_allMUArms1_TOI_end, {'-b', 'LineWidth', 1.3}, 1);
end

ylim([-0.5,0.5]);
if blAdj
title('Baseline Adjusted Average MUA RMS(1s) Aligned with Seizure Offset - Postictal and Recovery'); 
ylabel('MUA RMS(1s) Percentage Change');
else
title('Average MUA RMS(1s) Aligned with Seizure Offset');
ylabel('MUA RMS(1s)');
end
xlabel('Time (s)'); xlim([0,postIctalTime+recoveryTime]); 
xline(30,'-','Recovery','LabelHorizontalAlignment','right','LineWidth', 2);
set(MUA_RMS1_end, 'Position', [0,0,1900,600]);
fig_axs = findall(MUA_RMS1_end,'Type', 'Axes');
savenm = get(fig_axs(1), 'Title'); savenm = regexprep(savenm.String, ' ', '_'); 
saveas(MUA_RMS1_end, fullfile(svpath_ep, [savenm '.png']), 'png');
saveas(MUA_RMS1_end, fullfile(svpath_ep, [savenm '.eps']), 'eps');
saveas(MUA_RMS1_end, fullfile(svpath_ep, [savenm '.fig']), 'fig');

%% Stats
alpha = 0.05;

MUAallmean_mat = table2array(MUAallmean_table);
[pval, tbl, stats] = anova1(MUAallmean_mat);
if pval < alpha
    c = multcompare(stats, 'display', 'on', 'ctype','bonferroni');
    % finds a critical value for Dunnett's test by integrating the multivariate t distribution.
else
    % if anova does not return sig result, there is no point to run multcompare as no groups are significantly different from each other. 
    c = [];
end

MUAResults.Anova.p = pval;
MUAResults.Anova.table = tbl; 
MUAResults.Anova.stats = stats; 
MUAResults.Anova.MultComp = c; 

i = 1;
for pd_i = 1:4
    for pd_j = 1:4
        if pd_i < pd_j 
            % Signrank test
            signrankRes(i,1) = pd_i; signrankRes(i,2) = pd_j; 
            signrankRes(i,3) = signrank(MUAallmean_mat(:, pd_i), MUAallmean_mat(:, pd_j)); 
            % Paired t test
            ttestRes(i,1) = pd_i; ttestRes(i,2) = pd_j;
            [~,ttestRes(i,3)] = ttest(MUAallmean_mat(:, pd_i), MUAallmean_mat(:, pd_j)); 
            i = i +1;
        end
    end
end

MUAResults.Signrank = signrankRes;
MUAResults.PairedT = ttestRes;

%% save results
fprintf('Saving MUA analysis results ... \n');
save(fullfile(svpath_ep, 'MUAResAcrossPeriod.mat'),'MUAResults');
fprintf('Done!\n');


