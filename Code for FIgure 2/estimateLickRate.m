function lickRateTimeCourse = estimateLickRate(clickArray, afterLastClick, lickTimes, params, method)
% This function estimate the lick rate time course for each click from the
% lick events following that click. There are four ways to estimate this
% rate time course: (1) ordinary binning with a fixed bin size, (2) binning
% with variable bin size but fixed number of licks/bin, (3) sliding window
% function, (4) sliding gaussian kernel. 
% Inputs: 
%   - clickArray: array of click times for the clicks to be classified. 
%   - afterLastClick: time of the click after the last click so that we can
%   analyze the licks between the last click and click after it 
%   - lickTimes: array of all lick times in the seizure file to which the
%   click events belong. The time range of the clickArray must be contained
%   within the time range of the lickTimes array. This could also be a cell
%   array with the same number of entries as the number of clicks, in which
%   case each entry should contain the times of the licks for the
%   corresponding click. 
%   - params: a struct that contains at least the following fields:
%       > binSize
%       > preStimOffset
%       > postStimOffset
%       > SlidingWindowSize
%   - method: a string indicating which estimation method to use. could one
%   of the following:
%       > 'fixed_bin'
%       > 'variable_bin'
%       > 'sliding_window'
%       > 'sliding_kernel'
%
% Outputs: 
%   - lickRateTimeCourse: a cell where each entry corresponds to a click
%   and in each entry is a lick rate time course for that click. 
%
% Author: Abdo Sharaf - abdo.sharaf@yale.edu, Created: July 31 2020. Last Modified: July 31 2020. 
% *************************************************************************

%% get the lick events for each click if they're not provided directly (i.e. lickTimes argument is not a cell array) 
% lickTimes = Lick.times;
if ~iscell(lickTimes)
    % initialize the lickTimes cell array 
    lickTimesCell = cell(length(clickArray),1); 
    
    % loop through the click times and get the lick times for each click
    for clk = 1:length(clickArray)
        
        % get the index of the first lick and the total number of licks
        if clk < length(clickArray)
            % get the lick times for this click 
            lickTimesCell{clk} = lickTimes(lickTimes > clickArray(clk) & ...
                lickTimes < clickArray(clk+1), 1);
        else
            % get the lick times for this click 
            lickTimesCell{clk} = lickTimes(lickTimes > clickArray(clk) & ...
                lickTimes < afterLastClick, 1);
        end
        
        %%%%% xinyuan edit 27 April 2022
        if isempty(lickTimesCell{clk}) 
            lickTimesCell{clk} = 0;
        end
        %%%%%
        
    end
    
else
    lickTimesCell = lickTimes; 
end

%% estimate lick rate based on specified method 

% initialize the output array
lickRateTimeCourse = cell(length(clickArray),1); 

% get relevant parameters 
binSize = params.binSize; 
preStimOffset = params.PreStimOffset; 
postStimOffset = params.PostStimOffset; 

% array of time points indicating the start time (sec) of each bin (in the
% case of fixed_bin estimation)
binTimePts = -preStimOffset : binSize : postStimOffset - binSize; 

switch method
    case 'fixed_bin'
        % NOT YET IMPLEMENTED 
        
    case 'variable_bin'
        % NOT YET IMPLEMENTED 
        
    case 'sliding_window'
        
        interval = params.TimeCourseResolution;  
        windowSize = ceil(params.SlidingWindowSize/interval); 
        eventTimeCourses = cell(length(clickArray), 1); 
        
        for clk = 1:length(clickArray)
            
            clkLicks = lickTimesCell{clk}; 
            
            if clk < length(clickArray)
                timePts = clickArray(clk):interval:clickArray(clk+1);
            else
                timePts = clickArray(clk):interval:afterLastClick; 
            end
            
            tmcrs = zeros(length(timePts),1); 
            
            for tmpt = 1:length(tmcrs)
                tmcrs(tmpt) = numel(find(clkLicks> timePts(tmpt) & ...
                    clkLicks < timePts(tmpt)+interval)); 
            end
            
            eventTimeCourses{clk} = tmcrs; 
            
            % estimate rate with a sliding window 
            rateTmcrs = zeros(length(tmcrs), 1); 
            for pt = 1:length(tmcrs)
                
                startidx = uint8(pt - floor(windowSize / 2)); 
                endidx = uint8((pt + floor(windowSize / 2)) + uint8(mod(windowSize, 2) ~= 0)); 
                
                rateTmcrs(pt) = sum(tmcrs(max(startidx,1):min(endidx, length(tmcrs))))/(windowSize*interval); 
            end
            
            %%%%% xinyuan edit 27 April 2022
            lickRateTimeCourse{clk} = rateTmcrs; 
            if isempty(lickRateTimeCourse{clk})
                lickRateTimeCourse{clk} = 0;
            end
            %%%%%
                    
        end
end