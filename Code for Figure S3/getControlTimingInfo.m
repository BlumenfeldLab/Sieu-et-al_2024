function seizureInfo = getControlTimingInfo(seizurePath, seizureCell, params)
% This function reads a seizure's mat file and extracts information about
% the baseline times, control times, postcontrol times, recovery times. 
% It also extracts information about lick and click times. 
% INPUTS:
%   - seizurePath: path to the seizure file [string]
%   - seizureCell: cell array holding information about all seizures. [cell] 
%   - params: a struct that contains the following fields: [struct]
%               > SampleRate
%               > PreIctalOffset
%               > PostIctalOffset
%               > RecoveryDelay
%               > RecoveryPeriod 
%               > endTimeMethod
%               > electrodeLocation
%
% OUTPUTS:
%   - seizureInfo: a struct with the following fields:  [struct]
%               > AllTimes: an array with columns corresponding to time
%               stamps, each period has two columns for start and end 
%               > AllTimesIndcs: an array similar in structure to AllTimes
%               but is rounded appropriately to be used for indexing later.
%               > Lick: the lick information extracted from the file
%               > Click: the click information extracted from the file 
%               > ClickInMs: the click information in ms (samples)
%               > HC_LFP: the field with the relevant HC recording data 
%               > LO_LFP: the field with the relevant LO recording data
%               > hasControl: a logical flag to indicate whether the file
%               has a control period or not 
% Created by: Abdo Sharaf - abdo.sharaf@yale.edu
% Modified by: Xinyuan Zheng - xinyuan.zheng@yale.edu 
% Last updated: 02 Jul 2022 
%**************************************************************************


% seizurePath = filepath;
%% get the seizure filename and initialize parameters
% split the path into '\'-separated parts 
pathSplit = split(seizurePath,'\'); 

% get the file name
fileName = pathSplit{end}(1:end-4);
disp(fileName);
% get the seizure start and end times 
seizureInd = find(string(seizureCell(:,2))==fileName, 1);    % seizure index
startTime = seizureCell{seizureInd, 3}; 

switch params.endTimeMethod_beh
    case 'ipsi'
        endTime = seizureCell{seizureInd, 4}; 
    case 'contra'
        endTime = seizureCell{seizureInd, 5}; 
    case 'max_ipsi_contra'
        endTime = max([seizureCell{seizureInd,4}, seizureCell{seizureInd,5},seizureCell{seizureInd,8},seizureCell{seizureInd,9}]); 
    case 'add_delay'
        % add a delay to the ipsi end time
        tDiff = cell2mat(seizureCell(2:end,5)) - cell2mat(seizureCell(2:end,5)); 
        meanDiff = mean(tDiff);
        stdDiff = std(tDiff); 
        delay = meanDiff + 1 * stdDiff; 
        endTime = seizureCell(seizureInd, 4) + delay; 
    otherwise
        error([params.endTimeMethod_beh, ' is not a valid end time method.']); 
end

%% now define seizure times in ms
SeizureTimes = [startTime, endTime] * params.SampleRate; 
SeizureIndcs = round(SeizureTimes); 

%% load the seizure file and relevant fileds 
% load the mat file 
dataStruct = load(seizurePath); 
% get the names of the fields in the mat file 
ds_fieldNames = fieldnames(dataStruct); 
% get the index of the lick_on field 
lickOnIndx = contains(ds_fieldNames, 'lick_on'); 
% get the index of the dig mark field 
digMarkIndx = contains(ds_fieldNames, 'DigMark'); 
% get the index of the control time field 
ctrlTimeIndx = contains(ds_fieldNames, 'control_time'); 
% get the index of the LO LFP field 
loLFPInd = contains(ds_fieldNames, 'LO_LFP'); 
% get the index of the wheel field 
whlInd = contains(ds_fieldNames, 'wheel'); 
% get the index of the MUA VRMS field (if exists)
vrms0Ind = contains(ds_fieldNames, {'vrms0','VRMS0'}); 
% vrms0Ind = contains(ds_fieldNames, 'VRMS0'); 
vrms1Ind = contains(ds_fieldNames, {'vrms1','VRMS1'}); 
% vrms1Ind = contains(ds_fieldNames, 'VRMS1'); 
% get the raw mua data
muaInd = contains(ds_fieldNames, {'LO_MUA','MUA_LLO'}); 
% muaInd = contains(ds_fieldNames, 'MUA_LLO'); 

%% get the index of HC LFP recording based on the specified electrode location

% determine inclusion/exclusion based on the specified electrode
% location criteria
switch(params.electrodeLocation)
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
        possibleStrs = {'ipsi', 'ispi', 'ispsi', '000_HC_LFP',...
            '001_HC_LFP', '002_HC_LFP','ipLHC','ipRHC'};
    case 'contra'
        % possible strings corresponding to contra
        possibleStrs = {'contra','coRHC','coLHC'};
    case 'all'
        % possible strings corresponding to all
        possibleStrs = {'HC_LF', 'HC_ventral', 'HC_dorsal'}; % xinyuan edit
    otherwise
        error(['Invalid HC recording electrode location specification ''', hc_loc_select, ...
            '''. If you are specifying a new location, add it to the switch statement.']);
end
% find if a field contains the possible strings
electrodeLocationInd = contains(ds_fieldNames, possibleStrs);

%% Now load the data 
% convert the structure to a matrix just for ease of indexing 
dataCell = struct2cell(dataStruct); 
% get the licks 
Lick = dataCell{lickOnIndx}; 
% get the clicks 
Click = dataCell{digMarkIndx}; 
% get the clicks in Ms
ClickInMs = Click.times * params.SampleRate; 
% get the relevant HC recording field 
    % if there's more than one field, get the last one 
HC_LFP = dataCell{find(electrodeLocationInd, 1, 'last')}; 
% get the relevant LO recording field 
    % if there's more than one field, get the last one 
LO_LFP = dataCell{find(loLFPInd, 1, 'last')}; 
% get the wheel field 
hasWheel = 1; 
if any(whlInd)
    Wheel = dataCell{find(whlInd, 1)}; 
    seizureInfo.Wheel = Wheel; 
else
    fprintf(['Warning: The following seizure file does not have wheel data',...
        ': %s\n'], fileName); 
    hasWheel = 0; 
end
% mua vrms signal 
if any(vrms0Ind)
    vrms0 = dataCell{find(vrms0Ind, 1)}; 
    seizureInfo.VRMS0 = vrms0; 
end
if any(vrms1Ind)
    vrms1 = dataCell{find(vrms1Ind, 1)}; 
    seizureInfo.VRMS1 = vrms1; 
end

% mua raw signal
if any(muaInd)
    muaRaw = dataCell{find(muaInd, 1)}; 
    seizureInfo.MUA = muaRaw; 
end

%% define control
% If a control time exists, get it. Otherwise, or if it exists but is empty
% set assign it to Nan. 
hasControl = 1; 
if sum(ctrlTimeIndx)
    if ~isempty(dataCell{ctrlTimeIndx}.times)
        % get the seizure duration. the duration of the control period is
        % defined to be the same as the seizure duration
        seizureDur = params.SampleRate * (endTime - startTime); 
        % now define the control times 
        ControlTimes = [dataCell{ctrlTimeIndx}.times * params.SampleRate,...
            (dataCell{ctrlTimeIndx}.times * params.SampleRate)+seizureDur];
        ControlIndcs = [round(dataCell{ctrlTimeIndx}.times * params.SampleRate),...
            round(dataCell{ctrlTimeIndx}.times * params.SampleRate)+...
            round(seizureDur)]; 
    else
        ControlTimes = [NaN, NaN]; 
        ControlIndcs = ControlTimes; 
        hasControl = 0; 
    end
else
    ControlTimes = [NaN, NaN]; ControlIndcs = ControlTimes; 
    hasControl = 0; 
end



%% define baseline control 
baselineStart_ctrl = params.SampleRate * (dataCell{ctrlTimeIndx}.times - params.PreIctalOffset);
baselineStartInd_ctrl = round(ControlIndcs(1) - params.SampleRate*params.PreIctalOffset);
baseLineEnd_ctrl = params.SampleRate * dataCell{ctrlTimeIndx}.times - 1; 
baseLineEndInd_ctrl = round(ControlIndcs(1)) - 1; 

BaselineTimes_ctrl = [baselineStart_ctrl, baseLineEnd_ctrl];
BaselineIndcs_ctrl = [baselineStartInd_ctrl, baseLineEndInd_ctrl];

%% define postcontrol 

postStart_ctrl = (dataCell{ctrlTimeIndx}.times * params.SampleRate)+seizureDur + 1; 
postStartInd_ctrl = round((dataCell{ctrlTimeIndx}.times * params.SampleRate)+seizureDur) + 1; 

postEnd_ctrl = postStart_ctrl + params.SampleRate * params.PostIctalOffset - 1;
postEndInd_ctrl = postStartInd_ctrl + round(params.SampleRate * params.PostIctalOffset) - 1; 

PostTimes_ctrl = [postStart_ctrl, postEnd_ctrl]; 
PostIndcs_ctrl = [postStartInd_ctrl, postEndInd_ctrl]; 

%% define recovery control 

recoveryStart_ctrl = postEnd_ctrl + params.SampleRate * params.RecoveryDelay + 1; 
recoveryStartInd_ctrl = postEndInd_ctrl + round(params.SampleRate * params.RecoveryDelay) + 1; 

recoveryEnd_ctrl = recoveryStart_ctrl + params.SampleRate * params.RecoveryPeriod - 1; 
recoveryEndInd_ctrl = recoveryStartInd_ctrl + round(params.SampleRate * params.RecoveryPeriod) - 1; 

RecoveryTimes_ctrl = [recoveryStart_ctrl, recoveryEnd_ctrl]; 
RecoveryIndcs_ctrl = [recoveryStartInd_ctrl, recoveryEndInd_ctrl]; 

%% define the baseline period 
% start time of the baseline period
baseLineStart = params.SampleRate * (startTime - params.PreIctalOffset);
baseLineStartInd = round(startTime * params.SampleRate) - round(...
    params.PreIctalOffset * params.SampleRate); 
% end time of the baseline period 
baseLineEnd = params.SampleRate * startTime - 1; 
baseLineEndInd = round(startTime * params.SampleRate) - 1; 
% put them together 
BaselineTimes = [baseLineStart, baseLineEnd]; 
BaselineIndcs = [baseLineStartInd, baseLineEndInd]; 

%% define postictal times 
% start time of the postictal period
postIctalStart = params.SampleRate * endTime + 1; 
postIctalStartInd = round(params.SampleRate * endTime) + 1; 
% end time of the postictal period
postIctalEnd = postIctalStart + params.SampleRate * params.PostIctalOffset - 1;
postIctalEndInd = postIctalStartInd + round(params.SampleRate * params.PostIctalOffset) - 1; 
% put'em together
PostIctalTimes = [postIctalStart, postIctalEnd]; 
PostIctalIndcs = [postIctalStartInd, postIctalEndInd]; 

%% define recovery times 
% start time of the recovery period 
recoveryStart = postIctalEnd + params.SampleRate * params.RecoveryDelay + 1; 
recoveryStartInd = postIctalEndInd + round(params.SampleRate * params.RecoveryDelay) + 1; 
% end time of the recovery period 
recoveryEnd = recoveryStart + params.SampleRate * params.RecoveryPeriod - 1; 
recoveryEndInd = recoveryStartInd + round(params.SampleRate * params.RecoveryPeriod) - 1; 
% put em toghether
RecoveryTimes = [recoveryStart, recoveryEnd]; 
RecoveryIndcs = [recoveryStartInd, recoveryEndInd]; 



%% check if there's enough time for control

if ControlIndcs(2) < SeizureIndcs(1) % control end before seizure start
    
        % check if the start of the recording has enough time for baseline control 
    if baselineStart_ctrl < 0
        disp('err not enough time for baseline control(close to recording start)')
        hasControl = 0; 
    end
    
        % check if the post control is too close to the seizure baseline
    if recoveryEndInd_ctrl > baseLineStartInd
        disp('err not enough time for post control(close to seizure baseline start)')
        hasControl = 0; 
    end
    
elseif SeizureIndcs(2) < ControlIndcs(1) % seizure end before control start
    
    % check if the control start too close to recovery end, not enough time
    if baselineStart_ctrl < recoveryEndInd
        disp('err not enough time for baseline control(close to recovery end)')
        hasControl = 0; 
    end
        
    % check if recording is too short to have post control
    if recoveryEndInd_ctrl > length(LO_LFP.values)
        disp('err not enough time for post control(close to recording end)')
        hasControl = 0; 
    end
 
else
    disp('data err');
end


%% put everything into one array 
AllTimes = [BaselineTimes, SeizureTimes, PostIctalTimes, RecoveryTimes,...
    ControlTimes];
AllTimesIndcs = [BaselineIndcs, SeizureIndcs, PostIctalIndcs, RecoveryIndcs,...
    ControlIndcs]; 

AllTimes_ctrl = [BaselineTimes_ctrl, ControlTimes, PostTimes_ctrl, RecoveryTimes_ctrl];
AllTimesIndcs_ctrl = [BaselineIndcs_ctrl, ControlIndcs, PostIndcs_ctrl, RecoveryIndcs_ctrl]; 

%% put all relevant variables into the output struct 
seizureInfo.AllTimes = AllTimes; 
seizureInfo.AllTimesIndcs = AllTimesIndcs; 
seizureInfo.AllTimes_ctrl = AllTimes_ctrl;
seizureInfo.AllTimesIndcs_ctrl = AllTimesIndcs_ctrl; 
seizureInfo.Lick = Lick; 
seizureInfo.Click = Click; 
seizureInfo.ClickInMs = ClickInMs; 
seizureInfo.HC_LFP = HC_LFP; 
seizureInfo.LO_LFP = LO_LFP; 
seizureInfo.hasWheel = hasWheel; 
seizureInfo.hasControl = hasControl; 
end