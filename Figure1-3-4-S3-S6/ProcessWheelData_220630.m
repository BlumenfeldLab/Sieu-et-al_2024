% Edited by Jerry. 06/30/2022
function clickWheelData=ProcessWheelData_220630...
    (clickArray,wheelData,params)
    % ********************************************************************
    % This function sequences wheel data into windowed-timecourses
    % corresponding to individual clicks. 
    % It returns a matrix where rows are the number of samples in a window 
    % (specified in params by the preStimOffset and postStimOffset) and
    % the number of columns equals the number of clicks in the clickArray. 
    % Inputs:
    % (1) clickArray: array of click times 
    % (2) wheelData: the array with the wheel data values 
    % (3) params: the parameters struct from getParams function
    % Outputs: 
    % clickWheelData: a cell array with nClicks cells where in each cell
    %                 is the sequenced wheel data for each click
    % By Abdo Sharaf - abdo.sharaf@yale.edu
    % ********************************************************************
    %% Parameters initialization 
    % electrophysiology analysis pre and post stim times 
    preStim = params.PreStimOffset; 
    postStim = params.PostStimOffset; 
    numClicks = length(clickArray); % total number of clicks  
    fs = 1/wheelData.interval; % sample rate
    wheelData = wheelData.values; 
    % initialize output 
    clickWheelData = cell(numClicks, 1); 
    %% loop through clicks 
    for clk = 1:numClicks
       % get the start and end indices for the click's window 
       st = max(round((clickArray(clk)-preStim)*fs),1); 
       nd = min(round((clickArray(clk)+postStim)*fs),length(wheelData)); 
       % get the wheel data for this click 
       clickWheelData{clk} = wheelData(st:nd); 
    end
end