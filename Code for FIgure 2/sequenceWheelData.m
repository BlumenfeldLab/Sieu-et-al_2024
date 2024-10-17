function clickWheelData = sequenceWheelData(clickArray, wheelData, params)
% This function sequences wheel data into windowed-timecourses
% corresponding to individual clicks. It returns a matrix where the number
% of rows are the number of samples in a window (specified in params by the
% preStimOffset and postStimOffset) and the number of columns equals the
% number of clicks in the clickArray. 
% Inputs:
%   > clickArray: array of click times 
%   > wheelData: the array with the wheel data values 
%   > params: the parameters struct you get out of the getParams function
% Outputs: 
%   > clickWheelData: a cell array with nClicks cells where in each cell is 
%    the sequenced wheel data for each click
% By Abdo Sharaf - abdo.sharaf@yale.edu
% *************************************************************************

%% get relevant parameters 
% electrophysiology analysis pre and post stim times 
preStim = params.PreStimOffset; 
postStim = params.PostStimOffset; 
% total number of clicks 
numClicks = length(clickArray); 
% sample rate 
fs = 1/wheelData.interval; 
wheelData = wheelData.values; 
%% initialize output 
clickWheelData = cell(numClicks, 1); 

%% loop through clicks 

for clk = 1:numClicks
    
   % get the start and end indices for the click's window 
   st = max(round((clickArray(clk) - preStim)*fs), 1); 
   nd = min(round((clickArray(clk) + postStim)*fs), length(wheelData)); 
   
   % get the wheel data for this click 
   clickWheelData{clk} = wheelData(st:nd); 
   
end

end