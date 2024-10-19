% Edited by Jerry. 07/24/2022
%
% This function classifies click events as hits (spared seizures) or
% misses(impaired seizures) based on the classification method of choice.
% also returns the total number of licks associated with each click. 
%
function [clickClasses,featureArray,numLicks,timeCourses]=...
    AnalyzeClicks_220724...
    (clickArray,afterLastClick,lickTimes,binTimePts,params)
    %*********************************************************************
    % Inputs:
    % (1) clickArray: clicks' time points for clicks to be classified. 
    % (2) afterLastClick: the click's time after the last click so that we 
    %     can analyze lick respones between the last click and click after. 
    % (3) lickTimes: licks' times in the sz file to which clicks belong.
    %     clickArray's time range must be contained within lickTimes'
    %     time range, pass Lick.times to this function for classification
    % (4) binTimePts: the starting times (relative to a click) for bins
    % (5) params: a struct that contains at least the following fields:
    %     > firstLickWindow (classifying using the first lick method).
    %     > binSize
    %     > preStimOffset
    %     > postStimOffset
    %     > WinLickRateThreshold [classifying using the window lick rate
    %       method]
    %     > FractionOfBaseline [classifying using the normalized 
    %       window lick rate method)
    %
    % Outputs: 
    % (1) clickClasses: an array of the same length as the clickArray 
    %     where four possible values may exists:
    %     > (1) indicating the click event is a hit 
    %     > (0) indicating the click event is a miss
    %     > (-1) indicating the click event has no associated licks 
    %     > (2) indicating the click event is unclassified 
    % (2) numLicks: an array of the same length as the clickArray 
    %     that contains the lick number associated with each click.
    % (3) timeCourses: a matrix with nRows = nClicks and nCols = nBins, 
    %     which includes timecourses for all the clicks. 
    %     These are lick rate timecourses where the rate is estimated 
    %     using the binning procedure. 
    % (4) featureArray: features depending on the classification method.
    %     For example, feature array will be an array of response times
    %     if 'first_lick' classification method is used 
    %
    % CREATED BY: Abdo Sharaf - abdo.sharaf@yale.edu, ON: July 24th 2020
    %*********************************************************************
    %% parse optional arguments 
    % initialize click classes array, default as 'miss'
    clickClasses = zeros(length(clickArray),1);  
    % initialize the numLicks array 
    numLicks = zeros(length(clickArray),1); 
    % initialize the timecourses matrix 
    nBins = ceil((params.PreStimOffset+params.PostStimOffset)/...
        params.binSize); 
    timeCourses = zeros(length(clickArray),nBins); 
    % initialize the feature array
    featureArray = nan(length(clickArray), 1);
    %% analyze clicks based on 'first_lick'
    % loop through clicks 
    for click = 1:length(clickArray)
        % time of the click 
        clickTime = clickArray(click);  
        % get first lick index and the total number of licks
        if click < length(clickArray)
            % index of the first lick 
            firstLickInd = find(lickTimes>clickArray(click) & ...
                lickTimes<clickArray(click+1),1); 
            % number of licks 
            numLicks(click)=length(find(lickTimes>clickArray(click) & ...
                lickTimes<clickArray(click+1))); 
        else
            % index of the first lick 
            firstLickInd = find(lickTimes>clickArray(click) & ...
                lickTimes<=afterLastClick,1);
            % number of licks
            numLicks(click)=length(find(lickTimes>clickArray(click) & ...
                lickTimes<=afterLastClick)); 
        end
        % get the click's time course 
        for bin = 1:nBins
            % starting time corresponding the first bin edge
            st = clickTime + binTimePts(bin);   
            timeCourses(click,bin) = numel(find...
                ((lickTimes>st & lickTimes<(st+params.binSize)))); 
        end
        % if no associated licks
        if isempty(firstLickInd)            
            clickClasses(click) = -1;
        else
            % time of the first lick 
            firstLickTime = lickTimes(firstLickInd);  
            % delay to first lick 
            responseTime = firstLickTime - clickTime;  
            % classify based on firstLickWindow
            if responseTime <= params.firstLickWindow   
                clickClasses(click) = 1;  % hit (spared)
            end
            featureArray(click) = responseTime; 
        end
    end
end