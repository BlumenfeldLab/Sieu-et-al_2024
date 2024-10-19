function [MUAResponses, allClicks_mua, allClicks_lfp] = MUAPercentChange(seizureInfo, Behavior, params)
% this function extracts the vrms mua signal from the seizure files,
% divides it into clicks and normalizes it by baseline to compute percent
% change which can be later used to generate mua percent change plots for
% the hits and misses across different periods 
% By Abdo Sharaf - abdo.sharaf@yale.edu
% Xinyuan Zheng - xinyuan.zheng@yale.edu
% *************************************************************************

% seizureInfo = seizureInfos
%% Get relevant parameters 
% electrophysiology analysis pre and post stim times 
PreStimOffset = params.EPPreStimOffset; 
PostStimOffset = params.EPPostStimOffset; 
% seizure onset time skip 
szTimeSkip = params.SeizureOnsetTimeSkip;  
% total number of clicks 
numClicks = size(Behavior.ByClick.ClickTimes(:, 1), 1); 
% file numbers 
fileNums = sort(unique(Behavior.ByClick.ClickTimes(:, 1))); 
% sample rate 
SampleRate = params.SampleRate;
OtherFs = params.OtherSampleRate; 

%% initialize the variables that will hold analysis results 

% the click info 
ClickInfo.File = zeros(1, numClicks); 
ClickInfo.Class = nan(1, numClicks); 
ClickInfo.Skipped = nan(1, numClicks); 
ClickInfo.Period = zeros(1, numClicks);
ClickInfo.Onset = zeros(1, numClicks); 
ClickInfo.Animal = strings(1, numClicks);
ClickInfo.TimeFromPeriodOnset = zeros(1, numClicks); 

% mua responses
[MUAResponses.VRMS0, MUAResponses.VRMS1] = deal(nan((PreStimOffset+PostStimOffset)*SampleRate+1, numClicks));

% un-normalized sequenced signals
[allClicks_mua.vrms0, allClicks_mua.vrms1, allClicks_mua.raw] = deal(nan((PreStimOffset+PostStimOffset)*SampleRate+1, numClicks)); 
allClicks_lfp = nan((PreStimOffset+PostStimOffset)*OtherFs+1, numClicks); 

%% loop through the files in seizureInfo
clickCounter = 1; 
% array to hold the baseline averages for normalization
bslAvgs.vrms0 = zeros((PreStimOffset+PostStimOffset)*SampleRate+1, length(fileNums)); 
bslAvgs.vrms1 = zeros((PreStimOffset+PostStimOffset)*SampleRate+1, length(fileNums));

for file = 1:length(fileNums)
    fprintf('file number: %d \n', file);
    % get the relevant variables from inputs 
    %%% click times for the clicks in this file 
    fileClickInfo = Behavior.ByClick.ClickTimes(Behavior.ByClick.ClickTimes(:, 1)==fileNums(file),:); 
    fileClickTimes = fileClickInfo(:, 2); 
    
    %%% seizure onset and end in seconds 
    szStart = seizureInfo{file}.AllTimes(3)/SampleRate; 
    
    % get the vrms signal 
    vrms0 = seizureInfo{file}.VRMS0.values; 
    vrms1 = seizureInfo{file}.VRMS1.values; 
    
    % get the mua raw signal 
    mua = seizureInfo{file}.MUA.values; 
    
    % get the lfp signal 
    lfp = seizureInfo{file}.LO_LFP.values; 
    
    % skip clicks that are within the seizure onset skip period 
    skipCondition1 = fileClickTimes > szStart & (fileClickTimes - PreStimOffset) <= szStart + szTimeSkip; 
    % skip clicks that don't have enough pre stimulus data  &&
    skipCondition2 = round((fileClickTimes-PreStimOffset).*SampleRate) < 1;
    skipCondition3 = round((fileClickTimes+PostStimOffset).*SampleRate) > length(vrms0);
    
    clicksToSkip = skipCondition1 | skipCondition2 | skipCondition3;
    
    fileClks = 0;
    bslClks = 0;
    
    % loop through the clicks
    for clk = 1:length(fileClickTimes)
        % get the click's vrms signal
        if ~clicksToSkip(clk)
            
            % start and end indices 
            st = round(fileClickTimes(clk)*SampleRate) - round(PreStimOffset*SampleRate); 
            nd = round(fileClickTimes(clk)*SampleRate) + round(PostStimOffset*SampleRate); 
            stLFP = round(fileClickTimes(clk)*OtherFs) - round(PreStimOffset*OtherFs); 
            ndLFP = round(fileClickTimes(clk)*OtherFs) + round(PostStimOffset*OtherFs);
            % get the vrms0 for this click 
            %clkrms0 = vrms0(st:nd); 
            
            % get the vrms1 for this click 
            clkrms1 = vrms1(st:nd);
            
            % the raw signal
            clkmua = mua(st:nd); 
            
                        %%% TEMP %%% %% why??
            clkrms0 = movingrms(clkmua, 0.85, SampleRate);
            
            % get the lfp signal 
            clklfp = lfp(stLFP:ndLFP); 
            
            % append to the un-normalized matrix
            allClicks_mua.vrms0(:, clickCounter) = clkrms0; 
            allClicks_mua.vrms1(:, clickCounter) = clkrms1;
            allClicks_mua.raw(:, clickCounter) = clkmua; 
            allClicks_lfp(:, clickCounter) = clklfp; 
            
            fileClks = fileClks + 1; 
        end
        
        % populate the click information
        period = fileClickInfo(clk, 3);
        ClickInfo.File(clickCounter) = fileClickInfo(clk, 1); 
        ClickInfo.Class(clickCounter) = fileClickInfo(clk, 4);
        ClickInfo.Skipped(clickCounter) = clicksToSkip(clk);
        ClickInfo.Period(clickCounter) = period; 
        ClickInfo.Onset(clickCounter) = fileClickTimes(clk);
        ClickInfo.Animal(clickCounter) = Behavior.ByClick.AnimalsAll(clickCounter);
        
        if period == 1 || period == 2   % reference time for clicks in baseline is also the seizure start time
            ClickInfo.TimeFromPeriodOnset(clickCounter) = abs(fileClickTimes(clk) - szStart);
        else
            ClickInfo.TimeFromPeriodOnset(clickCounter) = fileClickTimes(clk) - seizureInfo{file}.AllTimes(2*period-1)/SampleRate;
        end
        
        % add to the baseline avg calculation of this file
        if period == 1
            bslAvgs.vrms0(:, file) = bslAvgs.vrms0(:, file) + clkrms0; 
            bslAvgs.vrms1(:, file) = bslAvgs.vrms1(:, file) + clkrms1; 
            bslClks = bslClks + 1;
        end
        
        clickCounter = clickCounter + 1; 
    end
    
    bslAvgs.vrms0(:, file) = bslAvgs.vrms0(:, file)./bslClks; 
    bslAvgs.vrms1(:, file) = bslAvgs.vrms1(:, file)./bslClks; 
    
    % normalize/calculate percent change wrt baseline averages 
    clkInds = ClickInfo.File == file; 
    MUAResponses.VRMS0(:, clkInds) = (allClicks_mua.vrms0(:, clkInds) - bslAvgs.vrms0(:, file))./bslAvgs.vrms0(:, file); 
    MUAResponses.VRMS1(:, clkInds) = (allClicks_mua.vrms1(:, clkInds) - bslAvgs.vrms1(:, file))./bslAvgs.vrms1(:, file);
end
MUAResponses.ClickInfo = ClickInfo; 
end