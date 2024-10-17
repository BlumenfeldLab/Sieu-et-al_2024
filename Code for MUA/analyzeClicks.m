function [clickClasses, featureArray, numLicks, timeCourses] = analyzeClicks(clickArray, afterLastClick,...
    lickTimes, binTimePts, params, varargin)
% this function classifies click events as hits (spared seizures) or
% misses(impaired seizures) based on the classification method of choice.
% also returns the total number of licks associated with each click. 
% Inputs:
%   - clickArray: array of click times for the clicks to be classified. 
%   - afterLastClick: time of the click after the last click so that we can
%   analyze the licks between the last click and click after it 
%   - lickTimes: array of all lick times in the seizure file to which the
%   click events belong. The time range of the clickArray must be contained
%   within the time range of the lickTimes array which is why it's best to
%   pass the Lick.times array directly to this function for classification 
%   - binTimePts: the starting times (relative to a click) for the bins
%   - wheelData: wheel data associated with the clicks. this will be used 
%   if we are classifiying clicks based on wheel data. [TODO]
%   which will be used to calculate the time courses 
%   - params: a struct that contains at least the following fields:
%       > firstLickWindow (in case of classifying using the delay to first
%       lick method). 
%       > binSize
%       > preStimOffset
%       > postStimOffset
%       > WinLickRateThreshold [in case of classifying using the window lick
%       rate method]
%       > FractionOfBaseline (in case of classifying using the normalized
%       window lick rate method)
%   - classificationMethod: this is an optional argument that indicates
%   which classification method to use. If not specified, the default
%   method is delay to first lick, for which a firstLickWindow parameter
%   must be passed. Two possible specifications for this argument are
%   supported as of now: 
%           >> 'first_lick': use the delay to first lick method 
%           >> 'window_lick_rate': use the window lick rate method 
%           >> 'win_lick_rate_and_first_lick': use both the delay to the
%           first lick and the window lick rate 
%           >> 'norm_lick_rate_and_first_lick': use both the delay to the
%           first lick and the baseline-normalized window lick rate. NOTE
%           that if this is specified, you have to pass in the baseline
%           click timcourse array (or in general the array you would like to
%           normalize by) as the last argument. the function will throw an
%           error if that's not the case 
%   - baselineTimeCourses: this is an optional argument that you have to pass if
%   you specified 'norm_lick_rate_and_first_lick' as your classification
%   method. this is the lick rate timecourses for the baseline clicks in 
%   by which you want to normalize the lick rate. 
% Outputs: 
%   - clickClasses: an array of the same length as the clickArray where
%   three possible values may exists:
%       > (1) indicating the click event is a hit 
%       > (0) indicating the click event is a miss
%       > (-1) indicating the click event has no associated licks 
%       > (2) indicating the click event is unclassified 
%   - numLicks: an array of the same length as the clickArray that contains
%   the number of licks associated with each click 
%   - timeCourses: a matrix with nRows = nClicks and nCols = nBins, which
%   includes timecourses for all the clicks. these are lick rate
%   timecourses where the rate is estimated using the binning procedure. 
%   - featureArray: an array of features depending on the classification
%   method used. For example, feature array will be an array of response
%   times in case 'first_lick' classification is used 
%
% CREATED BY: Abdo Sharaf - abdo.sharaf@yale.edu, ON: July 24th 2020
% EDITED BY: Xinyuan Zheng - xinyuan.zheng@yale.edu July 25 2022
%**************************************************************************

% lickTimes = Lick.times
%% parse optional arguments 
classMethod = 'first_lick'; 
if nargin > 3
    classMethod = varargin{1}; 
end

if strcmp(classMethod, 'norm_lick_rate_and_first_lick')
    if nargin < 7
        error(['Not enough input arguments for the classification method specified.',...
            ' Please, make sure you input the array by which to normalize']); 
    else
        normBy = varargin{2}; 
    end
end

% initialize click classes array 
clickClasses = zeros(length(clickArray), 1);    % they're all initialized as misses 

% initialize the numLicks array 
numLicks = zeros(length(clickArray), 1); 

% initialize the timecourses matrix 
nBins = ceil((params.PreStimOffset + params.PostStimOffset)/params.binSize); 
timeCourses = zeros(length(clickArray), nBins); 

%% analyze clicks  

switch classMethod
    case 'first_lick'
        % initialize the feature array
        featureArray = nan(length(clickArray), 1);
        
        % loop through the clicks 
        for click = 1:length(clickArray) % click = length(clickArray)
            
            clickTime = clickArray(click);  % time of the click 
            
            % get the index of the first lick and the total number of licks
            if click < length(clickArray)
                % index of the first lick 
                firstLickIndx = find(lickTimes > clickArray(click) & ...
                    lickTimes < clickArray(click+1), 1); 
                % number of licks 
                numLicks(click) = length(find(lickTimes > clickArray(click) & ...
                    lickTimes < clickArray(click+1))); 
            else
                firstLickIndx = find(lickTimes > clickArray(click) & ...
                    lickTimes <= afterLastClick, 1);
                numLicks(click) = length(find(lickTimes > clickArray(click) & ...
                    lickTimes <= afterLastClick)); 
            end
            
            % get the click's time course 
            for bin = 1:nBins
                st = clickTime + binTimePts(bin);   % starting time corresponding the first edge of the bin
                timeCourses(click, bin) = numel(find((lickTimes > st & lickTimes < (st + params.binSize)))); 
            end
            
            if isempty(firstLickIndx)           % if no associated licks 
                clickClasses(click) = -1;
                
            else
                firstLickTime = lickTimes(firstLickIndx);  % time of the first lick 
                responseTime = firstLickTime - clickTime;  % delay to first lick 
                
                % classify based on firstLickWindow
                if responseTime <= params.firstLickWindow   % hit (spared)
                    clickClasses(click) = 1; 
                end
                
                featureArray(click) = responseTime; 
            end
        end
        
    case 'window_lick_rate'  % this will classify based on the window lick rate (which is 
        % overall lick rate over a window of 0 to 2 seconds relative to
        % click onset)
        
        % initialize the feature array
        featureArray = nan(length(clickArray), 1);
        
        % loop through the clicks 
        for click = 1:length(clickArray)
            
            clickTime = clickArray(click);  % time of the click 
            
            % get the total number of licks
            if click < length(clickArray)
                
                % number of licks 
                numLicks(click) = length(find(lickTimes > clickArray(click) & ...
                    lickTimes < clickArray(click+1))); 
            else
                
                % number of licks
                numLicks(click) = length(find(lickTimes > clickArray(click) & ...
                    lickTimes <= afterLastClick)); 
            end
            
            % get the click's time course 
            for bin = 1:nBins
                st = clickTime + binTimePts(bin);   % starting time corresponding the first edge of the bin
                timeCourses(click, bin) = numel(find((lickTimes > st & lickTimes < (st + params.binSize)))); 
            end
            
            % calculate the window lick rate 
            stBinIdx = ceil(params.PreStimOffset / params.binSize);
            winLickRate = sum(timeCourses(click,stBinIdx:end), 2)/params.PostStimOffset;
            
            % add to the feature array 
            featureArray(click) = winLickRate; 
            
            % classify based on specified parameters 
            if winLickRate <= params.WinLickRateThreshold
                clickClasses(click) = 0;    % miss (impaired)
            else
                clickClasses(click) = 1;    % hit (spared)
            end
        end
        
    case 'win_lick_rate_and_first_lick'
        
        % initialize the feature array
        featureArray = nan(length(clickArray), 2);
        
        % loop through the clicks
        for click = 1:length(clickArray)
            
            clickTime = clickArray(click);  % time of the click
            
            % get the index of the first lick and the total number of licks
            if click < length(clickArray)
                % index of the first lick
                firstLickIndx = find(lickTimes > clickArray(click) & ...
                    lickTimes < clickArray(click+1), 1);
                
                % number of licks
                numLicks(click) = length(find(lickTimes > clickArray(click) & ...
                    lickTimes < clickArray(click+1)));
            else
                % index of the first lick
                firstLickIndx = find(lickTimes > clickArray(click) & ...
                    lickTimes <= afterLastClick, 1);
                
                % number of licks
                numLicks(click) = length(find(lickTimes > clickArray(click) & ...
                    lickTimes <= afterLastClick));
            end
            
            % get the click's time course
            for bin = 1:nBins
                st = clickTime + binTimePts(bin);   % starting time corresponding the first edge of the bin
                timeCourses(click, bin) = numel(find((lickTimes > st & lickTimes < (st + params.binSize))));
            end
            
            % calculate the window lick rate
            stBinIdx = ceil(params.PreStimOffset / params.binSize);
            winLickRate = sum(timeCourses(click,stBinIdx:end), 2)/params.PostStimOffset;
            
            % classify
            if isempty(firstLickIndx)           % if no associated licks
                % classify based on wind lick rate only
                if winLickRate <= params.WinLickRateThreshold
                    clickClasses(click) = 0;
                else
                    clickClasses(click) = 1;
                end
            else
                firstLickTime = lickTimes(firstLickIndx);  % time of the first lick
                responseTime = firstLickTime - clickTime;  % delay to first lick
                
                % classify based on both the delay and the window lick rate
                if responseTime <= params.firstLickWindow && winLickRate > params.WinLickRateThreshold      % hit (spared)
                    clickClasses(click) = 1;
                elseif responseTime > params.firstLickWindow && winLickRate <= params.WinLickRateThreshold    % miss (impaired)
                    clickClasses(click) = 0;
                else
                    clickClasses(click) = 2;    % 2 will correspond to unclassified
                end
                
                featureArray(click, 1) = responseTime;
            end
            
            featureArray(click, 2) = winLickRate;
        end
        
    case 'norm_lick_rate_and_first_lick'
        
        % initialize the feature array
        featureArray = nan(length(clickArray), 2);
        
        % get the average of the win lick rate for the baseline clicks by
        % which we are going to normalize 
        stBinIdx = ceil(params.PreStimOffset / params.binSize); % starting index for the window over which we calculate lick rate 
        winLickRateNorm = sum(normBy(:, stBinIdx:end), 2)./params.PostStimOffset; 
        normByValue = nanmean(winLickRateNorm); 

        
        % loop through the clicks
        for click = 1:length(clickArray)
            
            clickTime = clickArray(click);  % time of the click
            
            % get the index of the first lick and the total number of licks
            if click < length(clickArray)
                % index of the first lick
                firstLickIndx = find(lickTimes > clickArray(click) & ...
                    lickTimes < clickArray(click+1), 1);
                
                % number of licks
                numLicks(click) = length(find(lickTimes > clickArray(click) & ...
                    lickTimes < clickArray(click+1)));
            else
                % index of the first lick
                firstLickIndx = find(lickTimes > clickArray(click) & ...
                    lickTimes <= afterLastClick, 1);
                
                % number of licks
                numLicks(click) = length(find(lickTimes > clickArray(click) & ...
                    lickTimes <= afterLastClick));
            end
            
            % get the click's time course
            for bin = 1:nBins
                st = clickTime + binTimePts(bin);   % starting time corresponding the first edge of the bin
                timeCourses(click, bin) = numel(find((lickTimes > st & lickTimes < (st + params.binSize))));
            end
            
            % calculate the window lick rate
            winLickRate = sum(timeCourses(click,stBinIdx:end), 2)/params.PostStimOffset; 
            
            % normalize 
            normWinLickRate = winLickRate / normByValue; 
            
            % classify
            if isempty(firstLickIndx)           % if no associated licks
                % classify based on wind lick rate only
                if normWinLickRate <= params.FractionOfBaseline
                    clickClasses(click) = 0;
                else
                    clickClasses(click) = 1;
                end
            else
                firstLickTime = lickTimes(firstLickIndx);  % time of the first lick
                responseTime = firstLickTime - clickTime;  % delay to first lick
                
                % classify based on both the delay and the window lick rate
                if responseTime <= params.firstLickWindow && normWinLickRate > params.FractionOfBaseline      % hit (spared)
                    clickClasses(click) = 1;
                elseif responseTime > params.firstLickWindow && normWinLickRate <= params.FractionOfBaseline    % miss (impaired)
                    clickClasses(click) = 0;
                else
                    clickClasses(click) = 2;    % 2 will correspond to unclassified
                end
                
                featureArray(click, 1) = responseTime;
            end
            
            featureArray(click, 2) = normWinLickRate;
        end
      
    otherwise
        error('Invalid classification method sepcified'); 
end
                    
end