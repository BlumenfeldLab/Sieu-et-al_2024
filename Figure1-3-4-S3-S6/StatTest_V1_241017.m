% Fig 3E: LO delta Bandpower spared/impaired during ictal, postictal, late
% run DataAnalysis_V1_241017 with:
% AnimalID = 'Impairedsz',and 'Sparedsz';  HCstim = 'all'
% SavedFile\Impairedsz-all\DateTime\LODeltaBP.mat
% SavedFile\Sparedsz-all\DateTime\LODeltaBP.mat
%
% Fig 4E: LO delta Bandpower hit/miss during ictal
% run DataAnalysis_V1_241017 with:
% AnimalID = 'Allsz';  HCstim = 'all',
% SavedFile\Allsz-all\DateTime\ClkLO_DeltaBP.mat
% ClkLO_DeltaBP(2,:),ClkLO_DeltaBP(3,:)
%
% SupplFig 6C,D,G,H: HC delta/beta Bandpower spared/impaired during ictal
% run DataAnalysis_V1_241017 with:
% AnimalID = 'Impairedsz' and 'Sparedsz';  HCstim = 'ipsi' and 'contra'
% 
% SavedFile\Impairedsz-ipsi\DateTime\HCDeltaBP.mat
% SavedFile\Sparedsz-ipsi\DateTime\HCDeltaBP.mat
%
% SavedFile\Impairedsz-contra\DateTime\HCDeltaBP.mat
% SavedFile\Sparedsz-contra\DateTime\HCDeltaBP.mat
%
% SavedFile\Impairedsz-ipsi\DateTime\HCBetaBP.mat
% SavedFile\Sparedsz-ipsi\DateTime\HCBetaBP.mat
%
% SavedFile\Impairedsz-contra\DateTime\HCBetaBP.mat
% SavedFile\Sparedsz-contra\DateTime\HCBetaBP.mat

% normal distribution test using Shapiro-Wilk test or KS test
% if normally distributed, use unpaired ttest,
% if not, use Mann-Whitney U Test (Wilcoxon rank sum test)

% Edited by Jerry 10/17/2024
clear; 
close all;
clc;
%% dir for storing significance results 
% Initialization file saving directory
% SaveName = 'LODeltaBP-Sparedsz-Impairedsz-ictal';  % Fig 3E
% SaveName = 'LODeltaBP-Sparedsz-Impairedsz-postictal';  % Fig 3E
% SaveName = 'LODeltaBP-Sparedsz-Impairedsz-late';  % Fig 3E
% SaveName = 'ClkLO_DeltaBP-Hit-Miss-ictal'; % Fig 4E
% SaveName = 'ClkLO_DeltaBP-Hit-Miss-postictal'; % Fig 4E
SaveName = 'ClkLO_DeltaBP-Hit-Miss-late'; % Fig 4E
% SaveName = 'HCDeltaBP-Sparedsz-Impairedsz-ipsi';  % supplFig 6
% SaveName = 'HCDeltaBP-Sparedsz-Impairedsz-contra';  % supplFig 6
% SaveName = 'HCBetaBP-Sparedsz-Impairedsz-contra';  % supplFig 6
% SaveName = 'HCBetaBP-Sparedsz-Impairedsz-ipsi';  % supplFig 6
Datetime = datestr(now,'dd-mmm-yyyy_HH_MM_SS');
TempDir = ['E:\DataAnalysis\JerryLiu\SavedFile\' SaveName];
SaveDir = fullfile(TempDir,Datetime);
if ~exist(SaveDir,'dir')
    mkdir(SaveDir);
end
%% read data Fig 3E
% SaveTemp1 = 'E:\DataAnalysis\JerryLiu\SavedFile\Impairedsz-all\';
% SaveTemp2 = 'E:\DataAnalysis\JerryLiu\SavedFile\Sparedsz-all\';
% DataTemp1 = load([SaveTemp1 '29-Jul-2024_15_32_45\LODeltaBP.mat']);
% DataTemp2 = load([SaveTemp2 '29-Jul-2024_15_37_45\LODeltaBP.mat']);
% Data1 = DataTemp1.LODeltaBP;
% Data2 = DataTemp2.LODeltaBP;
%% read data Fig 4E
SaveTemp1 = 'E:\DataAnalysis\JerryLiu\SavedFile\JerryTest-all\';
SaveTemp2 = 'E:\DataAnalysis\JerryLiu\SavedFile\JerryTest-all\';
DataTemp1 = load([SaveTemp1 '29-Jul-2024_15_56_15\ClkLO_DeltaBP.mat']);
DataTemp2 = load([SaveTemp2 '29-Jul-2024_15_56_15\ClkLO_DeltaBP.mat']);
Data1 = DataTemp1.ClkLO_DeltaBP;
Data2 = DataTemp2.ClkLO_DeltaBP;
%% read data supplFig 6
% SaveTemp1 = 'E:\DataAnalysis\JerryLiu\SavedFile\Impairedsz-ipsi\';
% SaveTemp2 = 'E:\DataAnalysis\JerryLiu\SavedFile\Sparedsz-ipsi\';
% DataTemp1 = load([SaveTemp1 '29-Jul-2024_16_37_42\HCDeltaBP.mat']);
% DataTemp2 = load([SaveTemp2 '29-Jul-2024_16_44_13\HCDeltaBP.mat']);
% Data1 = DataTemp1.HCDeltaBP;
% Data2 = DataTemp2.HCDeltaBP;
%% read data supplFig 6
% SaveTemp1 = 'E:\DataAnalysis\JerryLiu\SavedFile\Impairedsz-contra\';
% SaveTemp2 = 'E:\DataAnalysis\JerryLiu\SavedFile\Sparedsz-contra\';
% DataTemp1 = load([SaveTemp1 '29-Jul-2024_16_42_22\HCDeltaBP.mat']);
% DataTemp2 = load([SaveTemp2 '29-Jul-2024_16_43_14\HCDeltaBP.mat']);
% Data1 = DataTemp1.HCDeltaBP;
% Data2 = DataTemp2.HCDeltaBP;
%% read data supplFig 6
% SaveTemp1 = 'E:\DataAnalysis\JerryLiu\SavedFile\Impairedsz-contra\';
% SaveTemp2 = 'E:\DataAnalysis\JerryLiu\SavedFile\Sparedsz-contra\';
% DataTemp1 = load([SaveTemp1 '29-Jul-2024_16_42_22\HCBetaBP.mat']);
% DataTemp2 = load([SaveTemp2 '29-Jul-2024_16_43_14\HCBetaBP.mat']);
% Data1 = DataTemp1.HCBetaBP;
% Data2 = DataTemp2.HCBetaBP;
%% read data supplFig 6
% SaveTemp1 = 'E:\DataAnalysis\JerryLiu\SavedFile\Impairedsz-ipsi\';
% SaveTemp2 = 'E:\DataAnalysis\JerryLiu\SavedFile\Sparedsz-ipsi\';
% DataTemp1 = load([SaveTemp1 '29-Jul-2024_16_37_42\HCBetaBP.mat']);
% DataTemp2 = load([SaveTemp2 '29-Jul-2024_16_44_13\HCBetaBP.mat']);
% Data1 = DataTemp1.HCBetaBP;
% Data2 = DataTemp2.HCBetaBP;
%% make sure Data1, Data2 are right vectors
% sData1 = Data1(2,:); % Fig 3E, supplFig ictal
% sData2 = Data2(2,:); % Fig 3E, supplFig ictal
% sData1 = Data1(3,:); % Fig 3E postictal
% sData2 = Data2(3,:); % Fig 3E postictal
% sData1 = Data1(4,:); % Fig 3E late
% sData2 = Data2(4,:); % Fig 3E late
% sData1 = Data1(2,:); % for click event, Fig 4E ictal
% sData2 = Data2(3,:); % for click event, Fig 4E ictal
% sData1 = Data1(4,:); % for click event, Fig 4E postictal
% sData2 = Data2(5,:); % for click event, Fig 4E postictal
sData1 = Data1(6,:); % for click event, Fig 4E late
sData2 = Data2(7,:); % for click event, Fig 4E late
%% normal distribution test using Shapiro-Wilk test or KS test
alpha = 0.01;
% [H1,pValue1,W1] = swtest_240729(swData1,alpha);
% [H2,pValue2,W2] = swtest_240729(swData2,alpha);
[H1,pValue1,W1] = kstest(sData1,'Alpha',alpha);
[H2,pValue2,W2] = kstest(sData2,'Alpha',alpha);
%% Unpaired t-test if normally distributed or Wilcoxon rank sum test. 
if H1&&H2 == 0
    [h,p,ci,stats] = ttest2(sData1,sData2,'Alpha',alpha,'Vartype','unequal');
    Ttest2Result.h = h; 
    Ttest2Result.p = p; 
    Ttest2Result.ci = ci;
    Ttest2Result.stats = stats;
    Title = ['Unpaired T-test-' num2str(alpha)];    
    save(fullfile(SaveDir,[Title '.mat']),'Ttest2Result'); 
else
    [p,h,stats] = ranksum(sData1,sData2,'Alpha',alpha);
    MWUtestResult.h = h; % 0: equal median, 1: unequal median
    MWUtestResult.p = p; % p>alpha: equal median, p<alpha: unequal median
    MWUtestResult.stats = stats; 
    Title = ['Mann-Whitney U-test-' num2str(alpha)]; 
    save(fullfile(SaveDir,[Title '.mat']),'MWUtestResult'); 
end
%% Nothing below