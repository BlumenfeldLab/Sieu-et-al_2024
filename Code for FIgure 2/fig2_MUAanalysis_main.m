%% Main script for the MUA analysis
% abdo.sharaf@yale.edu
% xinyuan.zheng@yale.edu 

%% paths 
clear; close all;

datapath = '\Data for Figure 2\Behavioral_partial_MUA'; 

svpath_beh = 'F:\xinyuan\Mouse Project\MUA Results\Behavior'; 
svpath_ep = 'F:\xinyuan\Mouse Project\MUA Results\EP';

t = datestr(now,'dd-mmm-yyyy_HH_MM_SS');
if ~exist(fullfile(svpath_beh,t), 'dir')
   mkdir(fullfile(svpath_beh,t));
end
if ~exist(fullfile(svpath_ep,t), 'dir')
   mkdir(fullfile(svpath_ep,t));
end

svpath_beh = fullfile(svpath_beh,t);
svpath_ep = fullfile(svpath_ep,t);

%% important variables
mua_fs = 20000; 
lfp_fs = 1000; 
names = string(ls(datapath)); 
names = names(3:end); 
saveFigs = true; 
params = getParams(); 
perAnimal = true; % can be either true or false
cls = {'Impaired', 'Spared'}; 
prds = {'Baseline', 'Ictal', 'PostIctal', 'Recovery'}; 

%% any desired changes to the params struct 
params.Analysis.MUAstats = 'average';    % can be changed to 'median'
params.electrodeLocation = 'all'; 
params.endTimeMethod = 'max_ipsi_contra'; 
params.ClassificationMethod = 'first_lick'; 
params.Normalization.Zero = 0; 
params.Analysis.ByAnimalFirst = false; 
params.Analysis.ByAnimalValue = 'average';
params.SampleRate = mua_fs; 
params.OtherSampleRate = lfp_fs; 
errstats = 'sem';

%% get seizure info and seizure paths 
fprintf('Getting information about all seizure files...\n'); 
[seizureCell, seizurePaths] = getAllSeizureInfo(names, datapath, params.electrodeLocation,false,svpath_beh);
fprintf('Done!\n'); 
animals = unique(seizureCell(2:end, 1)); 

%% behavior 
if ~perAnimal
    fprintf('Gathering behavioral information...\n'); 
    Behavior = LickTaskBehavior(seizurePaths, seizureCell, params, svpath_beh);
    fprintf('Done!\n'); 
else
    % loop through animals 
    Behavior = cell(length(animals), 1); 
    for an = 1:length(animals)
        animal = animals{an}; 
        animIndcs = find(strcmp(seizureCell, animal)); 
        miniCell = [seizureCell(1, :); seizureCell(animIndcs, :)];
        miniPaths = seizurePaths(animIndcs-1, :); 
        fprintf('Gathering behavioral information for %s ...\n', animal);
        if all(strcmp(miniCell(2,end-1), 'Exclude'))
            Behavior{an} = []; 
            fprintf('%s''s recordings have no clicks that satisfy selection parameters. Skipping...', animal); 
        else
            Behavior{an} = LickTaskBehavior(miniPaths, miniCell, params, svpath_beh); 
        end
        fprintf('Done!\n');         
    end
end

%% Looking at the raw data 
% loop through the files and get seizure infos
    for file = 1:length(seizurePaths)
    % path to the seizure file
    filepath = seizurePaths{file};
    % get the click, lick, and timing information for this seizure file
    szInfo = getSeizureTimingInfo(filepath, seizureCell, params);
    seizureInfos{file} = szInfo; 
    % get the seizure start time and indcs
    szStart = seizureCell{file+1, 3}; 
    startIdx_mua = ceil(szStart * mua_fs); 
    startIdx_lfp = ceil(szStart * lfp_fs);
    % get the lfp and mua signals 
    lfp = szInfo.LO_LFP.values(startIdx_lfp:end); 
    mua = szInfo.MUA.values(startIdx_mua:end); 
    % filter the signal 
    lfp_filtered = lowpass(lfp, 5, lfp_fs); 
    % downsample the mua signal and smooth the lfp signal 
    mua_ds = downsample(mua, mua_fs/lfp_fs); 
    lfp_smooth = smoothdata(lfp_filtered); 
    end

%% MUA Percent Change
if ~perAnimal
    fprintf('Gathering MUA Percent Change information...\n'); 
    [MUAResponses, allClicks_mua, allClicks_lfp] = MUAPercentChange(seizureInfos, Behavior, params);
    fprintf('Done!\n');     
else
   [MUAResponses, allClicks_mua, allClicks_lfp] = deal(cell(length(animals), 1)); 
   for an = 1:length(animals)
       animal = animals{an};
       animIndcs = find(strcmp(seizureCell, animal));
       miniSzInfos = seizureInfos(animIndcs-1); 
       fprintf('Gathering MUA Percent Change information for %s ...\n', animal);
       if ~isempty(Behavior{an})
           [MUAResponses{an}, allClicks_mua{an}, allClicks_lfp{an}] = MUAPercentChange(miniSzInfos, Behavior{an}, params);
       else
           [MUAResponses{an}, allClicks_mua{an}, allClicks_lfp{an}] = deal([]); 
       end
       fprintf('Done!\n');
   end
end


%% Plotting MUA Percent Change Data
figGroupNum = 0; figs = {}; 
if ~perAnimal
    figGroupNum = figGroupNum + 1; 
    for pd = 1:length(prds)
          % get all hits and misses
          clickinfo = MUAResponses.ClickInfo; 

          allhits = MUAResponses.VRMS0(:, clickinfo.Period == pd & clickinfo.Class == 1 & clickinfo.Skipped == 0);
          allmisses = MUAResponses.VRMS0(:, clickinfo.Period == pd & clickinfo.Class == 0 & clickinfo.Skipped == 0);

          % means and errors
          meanhits = nanmean(allhits, 2); 
          semhits = nanstd(allhits, [], 2)./sqrt(size(allhits, 2)); 
          meanmisses = nanmean(allmisses, 2); 
          semmisses = nanstd(allmisses, [], 2)./sqrt(size(allmisses, 2));
          
          % plot
          figs{figGroupNum}(pd) = figure('Visible', 'On');
          h1 = ShadedErrorBar((0:size(allhits, 1)-1)./mua_fs, meanhits,...
              semhits, {'-r', 'LineWidth', 1.3}, 1); hold on; 
          h2 = ShadedErrorBar((0:size(allmisses, 1)-1)./mua_fs, meanmisses,...
              semmisses, {'-b', 'LineWidth', 1.3}, 1);
          legend([h1.mainLine, h2.mainLine], ['Hits (n=', num2str(size(allhits, 2)), ')'], ['Misses (n=', num2str(size(allmisses, 2)), ')']);
          xlabel('Time (s)'); ylabel('VRMS'); 
          title(['Period = ', prds{pd}]); 
          
          % scatter plot of area under the curve
          allhits_auc = 0.001*trapz(allhits, 1);
          allmisses_auc = 0.001*trapz(allmisses, 1);
          if pd == 1
              data_cell{1} = allhits_auc'; 
              data_cell{2} = repmat({'Hits'}, length(allhits_auc'), 1); 
          else
              data_cell{1} = [allhits_auc'; allmisses_auc']; 
              data_cell{2} = [repmat({'Hits'}, length(allhits_auc'), 1); repmat({'Misses'}, length(allmisses_auc'), 1)]; 
          end
          figs{figGroupNum+1}(pd) = figure('Visible', 'On');
          boxplot(data_cell{1}, data_cell{2}); title(['AUC Across Hits and Misses for ' prds{pd} ' Period']); ylabel('AUC'); hold on
          scatter(gca, (1 + (rand(length(allhits_auc'), 1))), allhits_auc', 'r','filled', 'SizeData', 25, 'MarkerFaceColor',...
              'b', 'MarkerEdgeColor', 'none'); hold on
          scatter(gca, (2 + (rand(length(allmisses_auc'), 1))), allmisses_auc', 'r','filled', 'SizeData', 25, 'MarkerFaceColor',...
              'r', 'MarkerEdgeColor', 'none'); hold off
    end  
else
    % cell array to store the rms points 
    rms_animal = cell(1, length(animals)); 
    bsl_avgs = cell(1, length(animals)); 
    T = cell2table(cell(0,4));
    T.Properties.VariableNames = ["Animal","Period",'Nhits','Nmisses'];
    for an = 1:length(animals)
        figGroupNum = figGroupNum + 1;
        % get the baseline average of the rms 
        clickinfo = MUAResponses{an}.ClickInfo; 
        bsl_hits = allClicks_mua{an}.raw(:, clickinfo.Period == 1 & clickinfo.Class == 1 & clickinfo.Skipped == 0);
        avg_bslrms = nanmean(rms(bsl_hits, 1)); 
        bsl_avgs{an} = avg_bslrms; 
       for pd = 1:length(prds)
          % get all hits and misses
          allhits = MUAResponses{an}.VRMS0(:, clickinfo.Period == pd & clickinfo.Class == 1 & clickinfo.Skipped == 0);
          allmisses = MUAResponses{an}.VRMS0(:, clickinfo.Period == pd & clickinfo.Class == 0 & clickinfo.Skipped == 0);
          
          % get the raw hits and misses and their rms adjusted to baseline
          allhits_raw = allClicks_mua{an}.raw(:, clickinfo.Period == pd & clickinfo.Class == 1 & clickinfo.Skipped == 0);
          allmisses_raw = allClicks_mua{an}.raw(:, clickinfo.Period == pd & clickinfo.Class == 0 & clickinfo.Skipped == 0);
          avg_hitrms = nanmean(rms(allhits_raw, 1))/avg_bslrms -1; 
          avg_missrms = nanmean(rms(allmisses_raw, 1))/avg_bslrms -1; 
          rms_animal{an}(pd, 1) = avg_hitrms; 
          rms_animal{an}(pd, 2) = avg_missrms; 
          
          % means and errors
          meanhits = nanmean(allhits, 2); 
          semhits = nanstd(allhits, [], 2)./sqrt(size(allhits, 2)); 
          meanmisses = nanmean(allmisses, 2); 
          semmisses = nanstd(allmisses, [], 2)./sqrt(size(allmisses, 2));
          
          if 0 % plot
              figs{figGroupNum}(pd) = figure('Visible', 'On');
              h1 = ShadedErrorBar((0:size(allhits, 1)-1)./mua_fs, meanhits,...
                  semhits, {'-r', 'LineWidth', 1.3}, 1); hold on;
              h2 = ShadedErrorBar((0:size(allmisses, 1)-1)./mua_fs, meanmisses,...
                  semmisses, {'-b', 'LineWidth', 1.3}, 1);
              legend([h1.mainLine, h2.mainLine], ['Hits (n=', num2str(size(allhits, 2)), ')'], ['Misses (n=', num2str(size(allmisses, 2)), ')']);
              xlabel('Time (s)'); ylabel('VRMS');
              title(['Period = ', prds{pd}, ', Animal = ', animals{an}]);
              
              T = [T; {animals{an}, prds{pd}, num2str(size(allhits, 2)), num2str(size(allmisses, 2))}]
          end
       end
    end
end
writetable(T, ['Num_hits_misses_MUA_by_animal.xlsx'])

%% plot the average rms values across animals for each period (only if perAnimal is set to true)
if perAnimal
   figGroupNum = figGroupNum + 1; 
   for p = 2:length(prds)
       rms_array_balanced = nan(length(animals), 2); % this is with the constraint of an equal number of hits and misses
       rms_array_unbalanced = nan(length(animals), 2); % this is without that constraint 
       % loop through the animals and choose only the ones that have a
       % number for both hits and misses in this period
       for an = 1:length(animals)
          % exclude animals that don't have any clicks in at least one category (e.g. hits, misses)
          if ~isnan(rms_animal{an}(p, 1)) && ~isnan(rms_animal{an}(p, 2))
              rms_array_balanced(an, 1) = rms_animal{an}(p, 1); 
              rms_array_balanced(an, 2) = rms_animal{an}(p, 2);             
          end
          % include all animals 
          rms_array_unbalanced(an, 1) = rms_animal{an}(p, 1); 
          rms_array_unbalanced(an, 2) = rms_animal{an}(p, 2);
       end
       % convert to table 
       rms_table = array2table(rms_array_balanced, 'VariableNames', {'Hits', 'Misses'});
       
       if p == 2
           means = mean(rms_array_balanced,'omitnan');
           stds = std(rms_array_balanced,'omitnan');
           tmp = rms_array_balanced(:,1);
           sems = stds/sqrt(length(tmp(~isnan(tmp))));
           Fig4F_table = array2table([means;stds;sems], 'VariableNames', {'Hits', 'Misses'});
           writetable(Fig4F_table,'Fig4F_table.xlsx');
       end
       % now plot % the balanced plots
       figs{figGroupNum}(p-1) = figure('Visible', 'On'); 
       ttl1 = sprintf('Average Baseline-adjusted RMS of Each Animal for Hits and Misses in the %s Period (Balanced)', prds{p});
       boxplotScatter(rms_table, sprintf('Average Baseline-adjusted RMS of Each Animal for Hits and Misses in the %s Period (Balanced)', prds{p}),...
           'Baseline-adjusted Mean RMS', true, params.Analysis.MUAstats, errstats);
       set(gca,'fontsize',8);       set(gcf,'Position',[100 100 600 300]);
       saveas(gcf,fullfile(svpath_ep,[ttl1,'.png']))
       saveas(gcf,fullfile(svpath_ep,[ttl1,'.eps']))
       saveas(gcf,fullfile(svpath_ep,[ttl1,'.fig']))
       
       % the unbalanced plots 
       figs{figGroupNum+1}(p-1) = figure('Visible', 'On');  
       varNames = {'Hits', 'Misses'};
       if strcmp(params.Analysis.MUAstats,'median')
           boxplot(rms_array_unbalanced,'Labels',varNames,'symbol',''); hold on
       elseif strcmp(params.Analysis.MUAstats, 'average')
           if strcmp(errstats, 'std')
           errorbar([1:length(varNames)],mean(rms_array_unbalanced,'omitnan'),...
               std(rms_array_unbalanced,'omitnan'),'+','LineStyle','none',...
                'CapSize',18,'linewidth', 1,'MarkerSize',10,'color', 'k'); 
           elseif strcmp(errstats, 'sem')
               
           errorbar([1:length(varNames)],mean(rms_array_unbalanced,'omitnan'),...
               std(rms_array_unbalanced,'omitnan')./sqrt(length(sum(~isnan(rms_array_unbalanced),1))),...
               '+','LineStyle','none','CapSize',18,'linewidth', 1,'MarkerSize',10,'color', 'k'); 
           end
           xticks([1:length(varNames)]); xticklabels(varNames);hold on
       end
       
       colors = [10,139,148;130,10,40;134,179,82;242,110,48;177,224,123;232,214,16]./255;
       for hm = 1:2
           x = ones(size(rms_array_unbalanced(:,hm)))*hm;
           scatter(gca, x, rms_array_unbalanced(:,hm),colors(hm),'filled', 'SizeData', 25); hold on
       end
       
       hold off 
       xlim([0.8,2.2]);
       
       ttl2 = sprintf('Average Baseline-adjusted RMS of Each Animal for Hits and Misses in the %s Period (Unbalanced)', prds{p});
       title(sprintf('Average Baseline-adjusted RMS of Each Animal for Hits and Misses in the %s Period (Unbalanced)', prds{p})); 
       ylabel('Baseline-adjusted Mean RMS'); 
       set(gca,'fontsize',8);       set(gcf,'Position',[100 100 600 300]); 
       saveas(gcf,fullfile(svpath_ep,[ttl2,'.png']))
       saveas(gcf,fullfile(svpath_ep,[ttl2,'.eps']))
       saveas(gcf,fullfile(svpath_ep,[ttl2,'.fig']))

       % now do some stats 
       Results.RMS_Paired.(prds{p}).TwoTailedSignRank = signrank(rms_array_balanced(:, 1), rms_array_balanced(:, 2)); 
       Results.RMS_Paired.(prds{p}).RightTailedSignRank = signrank(rms_array_balanced(:, 1), rms_array_balanced(:, 2), 'tail', 'right');
       Results.RMS_Unpaired.(prds{p}).TwoTailedRankSum = ranksum(rms_array_unbalanced(:, 1), rms_array_unbalanced(:, 2));
       Results.RMS_Unpaired.(prds{p}).RightTailedRankSum = ranksum(rms_array_unbalanced(:, 1), rms_array_unbalanced(:, 2), 'tail', 'right');
   end
end

%% standardize and show plots
% standardizeScaleAndShow(figs, 1, svpath_ep);
close all;

%% save results
fprintf('Saving MUA analysis results ... \n');
save(fullfile(svpath_ep, 'MUAresults.mat'),'Results');
fprintf('Done!\n');
