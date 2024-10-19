% Loads animal recordings information. 
% Edited by Jerry, 09/20/2022
function [szCell,szDir,szTiming]=GetSeizureInfo_220919...
    (AnimalID,FileName,DataDir,params,bilateral,SaveDir) 
    % FUNCTION INFORMATION: Lick task analysis pipeline

    % Inputs:
    % (1)AnimalID: name of the animal [string].
    % (1)FileName: .mat data name [string].
    % (2)DataDir: .mat data directory [string].
    % (3)params: GetParameter output [struct]
    % (4)bilateral: true/false; true: exclude recordings which have no 
    %               bilateral implants. 
    % (5)SaveDir: Directory to which the szCell file be saved [string].

    % Outputs:
    % (1)szCell: cell with labelled seizure recording information.
    % (2)szDir: seizure-labelled .mat file directory
    %           (DataDir+AnimalID+Filename(szDataID)).
    % (3)szTiming: cell with baseline, sz, postictal, and recovery 
    %              start and end timing [s] and round[ms] index
    %              and events (lick, click, LFP, and wheel in structure)
    %              in four periods

    % szCell initialization
    seizureCellElements = {'animalID', 'szdataID','Include','Bilateral',...
        'startTime','endTime1','endTime2','endTime3','endTime4'};
    szCell = cell(size(seizureCellElements)); 
    szCell(1,:) = seizureCellElements; 
    % szDir initialization
    szDir = {};
    % szTiming initialization
    seizureTimingElements = {'AllTimes','AllTimesInds','Lick','Click',...
        'ClickInMs','HC_LFP','LO_LFP','hasWheel','Wheel','hasControl'};
    szTiming = cell(size(seizureTimingElements)); 
    szTiming(1,:) = seizureTimingElements; 
    % Go through each .mat file
    szCount = 0;
    for fn = 1:length(FileName(:,1))
        MatName = fullfile(DataDir,FileName(fn,:));
        disp(MatName);
        szDir{end+1,1}= MatName;
        % load seizure file 
        szInfo = load(MatName);
        % Bilateral recording?
        IsBilateral=sum(contains(fieldnames(szInfo),...
            {'ipsi','ispi','ispsi',...
            '000_HC_LF','001_HC_LF','002_HC_LF','ipLHC','ipRHC'}))&...
            sum(contains(fieldnames(szInfo),{'contra','coRHC','coLHC'}));
        % Inclusion/Exclusion based on electrode locations 
        switch(params.electrodeLocation)
            case 'ipsiLHC'
                % all possible strings corresponding to ipsiLHC
                possibleStrs = {'ipsiLHC_LFP','ipsiLHC_ventral_LFP',... 
                    'ipsiLHC_dorsal_LFP','ispiLHC_ventral_LFP',...  
                    'ispiLHC_dorsal_LFP','ispsiLHC_ventral_LFP',...  
                    'ispsiLHC_dorsal_LFP','ipLHC',...
                    '000_HC_LFP','001_HC_LFP','002_HC_LFP'};
            case 'ipsiRHC'
                % all possible strings corresponding to ipsiRHC
                possibleStrs = {'ipsiRHC_LFP','ipsiRHC_ventral_LFP',... 
                    'ipsiRHC_dorsal_LFP','ispiRHC_ventral_LFP',...  
                    'ispiRHC_dorsal_LFP','ispsiRHC_ventral_LFP',...  
                    'ispsRHC_dorsal_LFP','ipRHC'}; 
            case 'ipsi_all'
                % all possible strings corresponding to ipsi_all
                possibleStrs = {'ipsi','ispi','ispsi','ipRHC',... 
                    '001_HC_LF','002_HC_LF','ipLHC','000_HC_LF'}; 
            case 'contra'
                % possible strings corresponding to contra 
                possibleStrs = {'contra','coRHC','coLHC'};
            case 'all'
                % possible strings corresponding to all
                possibleStrs = {'HC_LF','HC_ventral','HC_dorsal'}; 
            otherwise
                error('Invalid HC electrode location specification.');
        end
        % find if a field contains the possible strings 
        electrodeLocationInd=contains(fieldnames(szInfo),possibleStrs);
        % determine exclusion and inclusion
        if ~IsBilateral&&bilateral
            toInclude = 'Exclude';
        % if a recording has an electrode in the desired location
        elseif sum(electrodeLocationInd)
            toInclude = 'Include'; 
        else
            toInclude = 'Exclude'; 
        end
        % get indices corresponding to seizure start and end fields
        StartInd = contains(fieldnames(szInfo),'seizure_time');
        EndInd1 = contains(fieldnames(szInfo),'seizure_end');
        EndInd2 = contains(fieldnames(szInfo),'seizure_2_end');
        % account for an update in the ipsilateral seizures end times
        EndInd3 = contains(fieldnames(szInfo),'seizure_3_end');  
        % account for an update in the contralateral seizures end times
        EndInd4 = contains(fieldnames(szInfo),'seizure_4_end');
        % get seizure start and end times 
        seizInfo = struct2cell(szInfo);
        szStart = seizInfo{StartInd}.times;
        szEnd1 = seizInfo{EndInd1}.times; % end time 1 ipsi endtime
        % if there is an end_2, get it as well 
        if sum(EndInd2)
            szEnd2 = seizInfo{EndInd2}.times; 
        else
            disp([string(MatName),'No endtime2'])
            szEnd2 = szEnd1;    % if not then assign it end1
        end
        % if there is an updated ipsi seizure time, update it as well
        if sum(EndInd3)
            szEnd3 = seizInfo{EndInd3}.times;
        else
            disp([string(MatName),'No endtime3'])
            szEnd3 = szEnd1;   
        end
        % if there is an updated ipsi seizure time, update it as well
        % szEnd3&szEnd4 is a pair during labelling process by Anna
        if sum(EndInd4)
            szEnd4 = seizInfo{EndInd4}.times;
        else
            disp([string(MatName),'No endtime4'])
            szEnd4 = szEnd3; 
        end
        % Transfer to string
        if IsBilateral
            Bilat = 'Yes'; 
        else
            Bilat = 'No'; 
        end
        % append info
        fnm = strtrim(FileName(fn,:)); 
        szCell(end+1,:)={AnimalID fnm(1:end-4) toInclude Bilat...
            szStart szEnd1 szEnd2 szEnd3 szEnd4};         
        % szTiming part below
        % check if this file is included or not based on Include
        % this means szTiming# = szCount and szCell# = fileCount
        if strcmp(toInclude,'Include')
            szCount = szCount + 1;  % increment the number of seizures 
        else
            continue;   % move on to next iteration 
        end
        % get endtime based on behavior or electrophysiology (EP)
        switch params.endTimeMethod
            case 'ipsi'
                szEnd = szEnd3; 
            case 'contra'
                szEnd = szEnd4; 
            case 'max_ipsi_contra'
                szEnd = max([szEnd1,szEnd2,szEnd3,szEnd4]); 
            otherwise
                error([params.endTimeMethod,...
                    ' end time method not valid.']);
        end
        % make it to ms 
        szTimes = [szStart,szEnd]*params.SampleRate; 
        szInds = round(szTimes);
        % get the names of the fields in the mat file 
        szInfoFieldNames = fieldnames(szInfo); 
        % get the index of the lick_on field 
        lickOnInd = contains(szInfoFieldNames,'lick_on'); 
        % get the index of the dig mark field 
        digMarkInd = contains(szInfoFieldNames,'DigMark'); 
        % get the index of the control time field 
        ctrlTimeInd = contains(szInfoFieldNames,'control_time'); 
        % get the index of the LO LFP field 
        lofLFPInd = contains(szInfoFieldNames,'LO_LFP'); 
        % get the index of the wheel field 
        whlInd = contains(szInfoFieldNames,'wheel'); 
        % get the index of the MUA VRMS field (if exists)
        vrms0Ind = contains(szInfoFieldNames,{'vrms0','VRMS0'}); 
        vrms1Ind = contains(szInfoFieldNames,{'vrms1','VRMS1'}); 
        % get the raw mua data
        muaInd = contains(szInfoFieldNames,{'LO_MUA','MUA_LLO'}); 
        % find if a field contains the possible strings
        electrodeLocationInd = contains(szInfoFieldNames,possibleStrs);
        % Convert structure to cell for easy index
        dataCell = struct2cell(szInfo); % cell elements are structures 
        Lick = dataCell{lickOnInd}; % Licks [structure]
        Click = dataCell{digMarkInd}; % Clicks [structure]
        ClickInMs = Click.times*params.SampleRate; % Clicks in ms [array]
        % Relevant HC recording fields,if more than one field,get the last 
        HC_LFP = dataCell{find(electrodeLocationInd,1,'last')}; 
        % Relevant LO recording field,if more than one field,get the last 
        LOF_LFP = dataCell{find(lofLFPInd,1,'last')}; 
        hasWheel = 1; % Wheel field 
        if any(whlInd)
            Wheel = dataCell{find(whlInd, 1)}; 
        else
            fprintf(['Warning: no wheel data in file',': %s\n'],...
                FileName(fn,:)); 
            hasWheel = 0; 
        end
        % mua vrms data 
        if any(vrms0Ind)
            vrms0 = dataCell{find(vrms0Ind, 1)}; 
        end
        if any(vrms1Ind)
            vrms1 = dataCell{find(vrms1Ind, 1)}; 
        end
        % mua raw data
        if any(muaInd)
            muaRaw = dataCell{find(muaInd, 1)}; 
        end
        % Get the control time, it no control time, assign it to Nan. 
        hasControl = 1; 
        if sum(ctrlTimeInd)
            if ~isempty(dataCell{ctrlTimeInd}.times)
                % get the seizure duration. the duration of the control 
                % period is defined to be the same as the seizure duration
                szDur = szEnd-szStart; %[s]
                % now define the control times 
                ControlTimes = [dataCell{ctrlTimeInd}.times*...
                    params.SampleRate,(dataCell{ctrlTimeInd}.times*...
                    params.SampleRate)+szDur*params.SampleRate];
                ControlInds = [round(dataCell{ctrlTimeInd}.times*...
                    params.SampleRate),round(dataCell{ctrlTimeInd}.times...
                    *params.SampleRate)+round(szDur*params.SampleRate)]; 
            else
                ControlTimes = [NaN,NaN]; 
                ControlInds = ControlTimes; 
                hasControl = 0; 
            end
        else
            ControlTimes = [NaN,NaN]; 
            ControlInds = ControlTimes; 
            hasControl = 0; 
        end
        % Define the baseline as the start of recroding until sz onset
        % Baseline period start
        baseLineStart = params.SampleRate*(szStart-params.PreIctalOffset);
        baseLineStartInd = round(szStart*params.SampleRate)-round(...
            params.PreIctalOffset*params.SampleRate); 
        % Baseline period end 
        baseLineEnd = params.SampleRate*szStart-1; 
        baseLineEndInd = round(szStart*params.SampleRate)-1; 
        % put them together 
        BaselineTimes = [baseLineStart, baseLineEnd]; 
        BaselineInds = [baseLineStartInd, baseLineEndInd];
        % Define postictal times 
        % Postictal period start
        postIctalStart = ...
            params.SampleRate*max([szEnd1,szEnd2,szEnd3,szEnd4])+1; 
        postIctalStartInd = ...
            round(params.SampleRate*max([szEnd1,szEnd2,szEnd3,szEnd4]))+1; 
        % Postictal period end
        postIctalEnd = postIctalStart+params.SampleRate*...
            params.PostIctalOffset-1;
        postIctalEndInd = postIctalStartInd+...
            round(params.SampleRate*params.PostIctalOffset)-1; 
        % put them together
        PostIctalTimes = [postIctalStart, postIctalEnd]; 
        PostIctalInds = [postIctalStartInd, postIctalEndInd]; 
        % Define recovery times 
        % Recovery period start
        recoveryStart = postIctalEnd+params.SampleRate*...
            params.RecoveryDelay+1; 
        recoveryStartInd = postIctalEndInd+...
            round(params.SampleRate*params.RecoveryDelay)+1; 
        % Recovery period end
        recoveryEnd = recoveryStart+params.SampleRate*...
            params.RecoveryPeriod-1; 
        recoveryEndInd = recoveryStartInd+...
            round(params.SampleRate*params.RecoveryPeriod)-1; 
        % put them toghether
        RecoveryTimes = [recoveryStart,recoveryEnd]; 
        RecoveryInds = [recoveryStartInd,recoveryEndInd];
        % put everything into one array 
        AllTimes = [BaselineTimes,szTimes,PostIctalTimes,RecoveryTimes,...
            ControlTimes];
        AllTimesInds = [BaselineInds,szInds,PostIctalInds,RecoveryInds,...
            ControlInds]; 
        % output
        szTiming(szCount,:)={AllTimes AllTimesInds Lick Click...
            ClickInMs HC_LFP LOF_LFP hasWheel Wheel hasControl};
%         szTiming(end+1,:)={AllTimes AllTimesInds Lick Click...
%             ClickInMs HC_LFP LOF_LFP hasWheel Wheel hasControl};
    end
    % Save szCell,szTiming
    save(fullfile(SaveDir,'seizureCell.mat'),'szCell'); 
%     save(fullfile(SaveDir,'seizureTiming.mat'),'szTiming'); 
    fprintf('Done!\n'); 
end