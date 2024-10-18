function [seizureCell, seizurePaths] = getAllSeizureInfo(subjectNames, dataPath, electrodeLocation, varargin)
% subjectNames = names;
% dataPath = datapath;
% electrodeLocation = params.electrodeLocation;
% bilateral = false;


% FUNCTION INFORMATION: Lick task analysis pipeline
% Created by: Abdo Sharaf - abdo.sharaf@yale.edu
% Xinyuan Zheng - xinyuan.zheng@yale.edu
% Date modified: 07/19/2022
%
% Loads information about animal recordings.
% Inputs:
% 1. subjectNames - string array of Animal folders names
% 2. dataPath - file path of the folder containing the animal folders/recording data.
% 3. electrodeLocation - string of HC electrode location to filter (for example,
%   if you wanted to get only ipsilateral (w.r.t side of stimulation
%   electrode) recordings that are on the left side, use 'ipsiLHC'. See the
%   switch statement for a list of valid inputs.
% 4. bilateral - if true, exclude recordings which do not have bilateral
%   implants. If you do not include this parameter, default value is false.
% 5. save_path - the path to which the seizureCell file should be saved. If
% not specified, the seizureCell file will be saved to the provided data
% path by default. 
%
% Outputs:
% 1. seizureCell - cell with recording information. The features included are
%   described in the cell 'seizureCellElements' below. This cell is also
%   saved as a .mat file.
%
% 2. seizurePaths - string array with full paths to seizures
%*********************************************************************************

% Set default value of bilateral to false if input not provided
if nargin < 4
    bilateral = false; 
    savePath = dataPath; 
elseif nargin < 5
    bilateral = varargin{1}; 
    savePath = dataPath;
elseif nargin >= 5
    bilateral = varargin{1}; 
    savePath = varargin{2}; 
end

% initialize the cell with headings (also used as descriptors for features
% to select
% seizureCellElements = {'animalID', 'seizureID', 'startTime', 'endTime',...
%     'endTime2', 'Include', 'Bilateral'};
seizureCellElements = {'animalID', 'seizureID', 'startTime', 'endTime',...
    'endTime2','Include', 'Bilateral','endTime3','endTime4','HasControl'};
seizureCell = cell(size(seizureCellElements)); 
seizureCell(1,:) = seizureCellElements; 

% initialize the cell array with the seizure paths 
seizurePaths = {}; 

fprintf('Getting information about all seizure files...\n'); 
% Loop through all recording files and get appropriate variables 
for fl = 1:length(subjectNames)
    % path to the folder of the current animal 
    currSubjectFolder = fullfile(dataPath, subjectNames(fl)); 
    % get names of all the seizure files within this animal's folder 
    seizFiles = ls(fullfile(currSubjectFolder,'*.mat'));
    
    % loop through all seizure files and extract relevant info 
    for sz = 1:size(seizFiles, 1)
        % sz = 4
        
        % path to the seizure file 
        seizPath = fullfile(currSubjectFolder, string(seizFiles(sz, :)));
        seizurePaths{end+1, 1} = seizPath;  % append to the paths cell
        
        disp(seizPath);
        
        % load seizure file 
        seizInfo = load(seizPath); 
        
        % is this a bilateral recording? % xinyuan edited 4 April 2022
        is_bilateral = sum(contains(fieldnames(seizInfo),{'ipsi', 'ispi',...
            'ispsi', '000_HC_LF','001_HC_LF', '002_HC_LF','ipLHC','ipRHC'})) & ...
            sum(contains(fieldnames(seizInfo), {'contra','coRHC','coLHC'})); % TODO 
        
        
        % determine inclusion/exclusion based on the specified electrode
        % location criteria 
        switch(electrodeLocation)
            case 'ipsiLHC'
                % all possible strings corresponding to ipsiLHC
                possibleStrs = {'ipsiLHC_LFP', 'ipsiLHC_ventral_LFP',... 
                    'ipsiLHC_dorsal_LFP', 'ispiLHC_ventral_LFP',...  
                    'ispiLHC_dorsal_LFP', 'ispsiLHC_ventral_LFP',...  
                    'ispsiLHC_dorsal_LFP', '000_HC_LFP', '001_HC_LFP',... 
                    '002_HC_LFP','ipLHC'};  
            case 'ipsiRHC'
                % all possible strings corresponding to ipsiRHC
                possibleStrs = {'ipsiRHC_LFP', 'ipsiRHC_ventral_LFP',... 
                    'ipsiRHC_dorsal_LFP', 'ispiRHC_ventral_LFP',...  
                    'ispiRHC_dorsal_LFP', 'ispsiRHC_ventral_LFP',...  
                    'ispsRHC_dorsal_LFP','ipRHC'}; 
            case 'ipsi_all'
                % all possible strings corresponding to ipsi_all
                possibleStrs = {'ipsi', 'ispi', 'ispsi', '000_HC_LF',... 
                    '001_HC_LF', '002_HC_LF','ipLHC','ipRHC'}; 
            case 'contra'
                % possible strings corresponding to contra 
                possibleStrs = {'contra','coRHC','coLHC'};
            case 'all'
                % possible strings corresponding to all
                possibleStrs = {'HC_LF', 'HC_ventral', 'HC_dorsal'};  % xinyuan edit for FP data 4 April 2022
            otherwise
                error(['Invalid HC recording electrode location specification ''', hc_loc_select, ...
                    '''. If you are specifying a new location, add it to the switch statement.']);
        end
        % find if a field contains the possible strings 
        electrodeLocationInd = contains(fieldnames(seizInfo), possibleStrs);
        
        % determine exclusion and inclusion
        if ~is_bilateral && bilateral
            toInclude = 'Exclude'; 
        elseif sum(electrodeLocationInd)    % if a recording has an electrode in the desired location
            toInclude = 'Include'; 
        else
            toInclude = 'Exclude'; 
        end
        
        % get indices corresponding to seizure start and end fields
        start_ind = contains(fieldnames(seizInfo), 'seizure_time');
        end_ind = contains(fieldnames(seizInfo), 'seizure_end');
        end_ind2 = contains(fieldnames(seizInfo), 'seizure_2_end');
        % account for an update in the ipsilateral seizures end times  % xinyuan 02 May 2022
        end_ind3 = contains(fieldnames(seizInfo), 'seizure_3_end');  
        % account for an update in the contralateral seizures end times
        end_ind4 = contains(fieldnames(seizInfo), 'seizure_4_end');
        
        ctrlTimeIndx = contains(fieldnames(seizInfo), 'control_time'); 

        % get seizure start and end times 
        seizInfo = struct2cell(seizInfo);
        
        % If a control time exists, get it. Otherwise, or if it exists but is empty
        % set assign it to Nan. 
        hasControl = 1; 
        if sum(ctrlTimeIndx)
            if isempty(seizInfo{ctrlTimeIndx}.times)
                hasControl = 0; 
            end
        else
            hasControl = 0; 
        end
        
        seizStart = seizInfo{start_ind}.times;
        
        % end time 1 ipsi endtime
        seizEnd1 = seizInfo{end_ind}.times; 
        
        % if there is an end_2 get it as well 
        if sum(end_ind2)
            seizEnd2 = seizInfo{end_ind2}.times; 
        else
            disp([string(seizFiles(sz, :)),'no endtime2 '])
            seizEnd2 = seizEnd1;    % if not then assign it the same value as end1
        end
        
        % if there is an updated ipsi seizure time, update it here as well
        if sum(end_ind3)
            seizEnd3 = seizInfo{end_ind3}.times;
        else
            disp([string(seizFiles(sz, :)),'no seizEnd3 '])
            seizEnd3 = seizEnd1;   
        end
        
        % if there is an updated ipsi seizure time, update it here as well
        if sum(end_ind4)
            seizEnd4 = seizInfo{end_ind4}.times;
        else
            disp([string(seizFiles(sz, :)),'no seizEnd4 '])
            seizEnd4 = seizEnd3;
        end
        
        if is_bilateral
            bilat = 'Yes'; 
        else
            bilat = 'No'; 
        end
        
        % append info
        fnm = strtrim(seizFiles(sz,:)); 
%         seizureCell(end+1,:) = {char(strtrim(subjectNames(fl))) fnm(1:end-4)...
%             seizStart seizEnd3 seizEnd4 toInclude bilat}; 
        seizureCell(end+1,:) = {char(strtrim(subjectNames(fl))) fnm(1:end-4)...
            seizStart seizEnd1 seizEnd2 toInclude bilat seizEnd3 seizEnd4 hasControl}; % try different ends / plot xinyuan 4 May 2022
    end
end

%% Save 
save(fullfile(savePath,'seizureCell.mat'), 'seizureCell'); 
fprintf('Done!\n'); 
end