% Main script for the Behavior analysis for Control periods
% By: xinyuan Zheng Jan 18 2022 - xinyuan.zheng@yale.edu 
%% paths 
clear; close all;

group = 'BehControl_1uA';
datapath = 'F:\xinyuan\Mouse Project\Behavior_partial_confirmed_LOVO'; 
res_path = 'F:\xinyuan\Mouse Project\Results'; 

%% important variables 
names = string(ls(datapath)); 
names = names(3:end); 
saveFigs = true; 
params = getParams(); 
perAnimal = false; 

%% any desired changes to the params struct 
params.electrodeLocation = 'all'; 
params.endTimeMethod_beh = 'max_ipsi_contra'; % use 'ipsi': end1, use 'contra': end2
params.ClassificationMethod = 'first_lick'; 
params.Normalization.Zero = 0; 
params.Analysis.ByAnimalFirst = false; 
params.Analysis.ByAnimalValue = 'average';

%%
t = datestr(now,'dd-mmm-yyyy_HH_MM_SS'); 
if ~exist(fullfile(res_path,[t,'_',group,'_Beh_end',params.endTimeMethod_beh]), 'dir')
   mkdir(fullfile(res_path,[t,'_',group,'_Beh_end',params.endTimeMethod_beh]));
end
svpath_beh = fullfile(res_path,[t,'_',group,'_Beh_end',params.endTimeMethod_beh]);

%% get seizure info and seizure paths 
[seizureCell, seizurePaths] = getAllSeizureInfo(names, datapath, params.electrodeLocation, false, svpath_beh);
seizureCell_ctrl = seizureCell([1; find(cell2mat(seizureCell(2:end,10)))+1],:);
seizurePaths_ctrl = seizurePaths(find(cell2mat(seizureCell(2:end,10))));
% check if has enough time for control analysis
toInclude = string(seizureCell_ctrl(2:end,6));
for file = 1:length(seizurePaths_ctrl) % file = 116
    filepath = seizurePaths_ctrl{file};
    disp(filepath);
    % check if this file should be included or not. 
    if ~strcmp(toInclude(file), 'Include')
        continue;   % move on to next iteration
    end
    % get the click, lick, and timing information for this seizure file
    seizureInfo = getControlTimingInfo(filepath, seizureCell_ctrl, params);
    seizureCell_ctrl{file+1,end} = seizureInfo.hasControl;
end
seizureCell_ctrl_rm = seizureCell_ctrl([1; find(cell2mat(seizureCell_ctrl(2:end,10)))+1],:);
seizurePaths_ctrl_rm = seizurePaths_ctrl(find(cell2mat(seizureCell_ctrl(2:end,10))));

%% behavior 
seizureCell = seizureCell_ctrl_rm;
seizurePaths = seizurePaths_ctrl_rm;
Behavior = LickTaskBehavior_ctrl(seizurePaths,seizureCell, params, svpath_beh);
%% analyzing behavior (not yet implemented for perAnimal analysis)
behResults = analyzeBehavior_ctrl(Behavior, params, svpath_beh); 
save(fullfile(svpath_beh, 'behResults.mat'),'behResults');
close all;
disp('Done!')

