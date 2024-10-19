% Edited by Jerry. 06/30/2022
%
% This function estimate the lick rate time course for each click from the
% lick events following that click. There are four ways for this estimate:
% (1) ordinary binning with a fixed bin size,
% (2) binning with variable bin size but fixed number of licks/bin,
% (3) sliding window function,
% (4) sliding gaussian kernel. 
% The output of this function did not be used
function lickRateTimeCourse = EstimateLickRate_220630...
    (clickArray,afterLastClick,lickTimes,params)
    % ********************************************************************
    % Inputs: 
    % (1) clickArray: clicks' time points for clicks to be classified. 
    % (2) afterLastClick: the click's time after the last click within
    %     a period so we can analyze lick respones between the last click
    %     in the period and the first click in the next period.                
    % (3) lickTimes: licks' times in the sz file to which clicks belong.
    %     clickArray's time range must be contained within lickTimes'
    %     time range, pass Lick.times to this function for classification
    % (4) params: a struct that contains at least the following fields:
    %     > binSize
    %     > preStimOffset
    %     > postStimOffset
    %     > SlidingWindowSize
    % ********************************************************************
    % Outputs: 
    %   lickRateTimeCourse: a cell where each entry corresponds to a click
    %   and in each entry is a lick rate time course for that click. 
    % ********************************************************************
    % Author: Abdo Sharaf - abdo.sharaf@yale.edu, ...
    % Created: July 31 2020. Last Modified: July 31 2020. 
    % ********************************************************************
    %% get the lick events for each click if they're not provided directly 
    % (i.e. lickTimes argument is not a cell array) 
    if ~iscell(lickTimes)
        % initialize the lickTimes cell array 
        lickTimesCell = cell(length(clickArray),1); 
        % loop through click times and get lick times for each click
        for clk = 1:length(clickArray)
            % get the index of the first lick and the total number of licks
            if clk < length(clickArray)
                % get the lick times for this click 
                lickTimesCell{clk}=lickTimes(lickTimes>clickArray(clk)&...
                    lickTimes<clickArray(clk+1), 1);
            else
                % get the lick times for this click 
                lickTimesCell{clk}=lickTimes(lickTimes>clickArray(clk)&...
                    lickTimes<afterLastClick, 1);
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
    lickRateTimeCourse = cell(length(clickArray),1); % output initialize
    interval = params.TimeCourseResolution;  % 0.01s
    windowSize = ceil(params.SlidingWindowSize/interval); % 0.1/0.01=10
    eventTimeCourses = cell(length(clickArray), 1); 
    for clk = 1:length(clickArray)
        clkLicks = lickTimesCell{clk}; 
        if clk < length(clickArray)
            timePts = clickArray(clk):interval:clickArray(clk+1);
        else
            timePts = clickArray(clk):interval:afterLastClick; 
        end
        %
        tmcrs = zeros(length(timePts),1); % # of licks within each interval
        for tmpt = 1:length(tmcrs)
            tmcrs(tmpt) = numel(find(clkLicks>timePts(tmpt) & ...
                clkLicks<timePts(tmpt)+interval)); 
        end
        eventTimeCourses{clk} = tmcrs; 
        % estimate rate with a moving window with step = 0.01s 
        rateTmcrs = zeros(length(tmcrs),1); 
        for pt = 1:length(tmcrs)
            startidx = uint8(pt-floor(windowSize/2)); 
            endidx = uint8((pt+floor(windowSize/2))+...
                uint8(mod(windowSize,2)~=0)); 
            % lick rate distribution with 0.1s resolution
            rateTmcrs(pt) = sum...
                (tmcrs(max(startidx,1):min(endidx,length(tmcrs))))/...
                (windowSize*interval); 
        end
        %%%%% xinyuan edit 27 April 2022
        lickRateTimeCourse{clk} = rateTmcrs; % 0.1s resolution lick rate
        if isempty(lickRateTimeCourse{clk})
            lickRateTimeCourse{clk} = 0;
        end
        %%%%%
    end
end