function params = getParams(varargin)
% this function get the suitable parameters for analysis. call this
% function before implementing the analysis and pass the output struct as
% an argument to functions that need it. if you are to change any of the
% parameters, change them here, not within functions. The only optional
% argument here is the save path to which to save the params mat file. If a
% save path is not specified, will save to the MATLAB path. 
% Here are the parameters defined here:
%   > SampleRate
%   > PreIctalOffset
%   > PostIctalOffset
%   > RecoveryDelay
%   > RecoveryPeriod 
%   > endTimeMethod
%   > electrodeLocation
%   > firstLickWindow
%   > maxLickRateWindow
%   > binSize
%   > PreStimOffset
%   > PostStimOffset
%   > trialWindow
%   > TimeCourseResolution 
%   > SlidingWindowSize
%   > WinLickRateThreshold
%   > ClassificationMethod
%   > SpectrogramWindow
%   > SpectrogramOverlap
%   > StimulationArtifactPeriod
%   > SpectrogramFrequencyRange
% Created by Abdo Sharaf - abdo.sharaf@yale.edu, July 25th 2020
% Xinyuan Zheng - edit June 2 2022
% *************************************************************************

%% parse optional args 
savePath = pwd; 
if nargin > 0
    savePath = varargin{1}; 
end

%% define the parameters    
params.SampleRate = 1000;                        % the sampling rate 

params.BaselineRef = true;                       % whether to regard all clicks in baseline as hits (or reference clicks really)
        
params.PreIctalOffset = 60;                      % duration in seconds of the pre ictal (baseline) period 

params.PostIctalOffset = 30;                     % duration in seconds of the post-ictal period 

params.RecoveryDelay = 0;                        % duration in seconds of the delay between the post-ictal period and the recovery period 

params.RecoveryPeriod = 60;                      % duration in seconds of the recovery period 

params.endTimeMethod_beh = 'max_ipsi_contra';    % if looking at both ipsi and contra, which end time should we use for Behavior analysis

params.endTimeMethod_EP = 'ipsi';                % which end time should we use for EP analysis

params.electrodeLocation = 'all';                % electrode location of interest, other locations are discarded 

params.firstLickWindow = 1.05;                   % threshold for the delay to first lick used for when classifying with first_lick method 

params.maxLickRateWindow = [-0.1, 1.05];         % interval (time before and time after the click) within which to calculate max lick rate 

params.binSize = 0.1;                            % the size of the bin used for estimating lick rate with binning 

params.PreStimOffset = 2;                        % for how long before the click should we analyze behavior
    
params.PostStimOffset = 2+0.1;                   % for how long after the stim should we analyze behavior

params.EPPreStimOffset = 2;                      % for how long before the click should we analyze electrophysiology 

params.EPPostStimOffset = 2;                     % for how long after the click should we analyze electrophysiology 

params.trialWindow = 3;                          % time window for the click beyond which responses are considered timed out     

params.TimeCourseResolution = 0.01;              % interval between consecutive time points in the time course of estimated lick rate 

params.SlidingWindowSize = 0.1;                  % size of the sliding window in when estimating the lick rate with a sliding window method 

params.WinLickRateThreshold = 0.5;               % threshold used for when classifiying based on window lick rate 

params.ClassificationMethod = 'first_lick';      % click classification method 

params.FractionOfBaseline = 0.5;                 % fraction of baseline to use as criteria when classifying with normalized lick rate

params.SpectrogramWindow = 1000;                 % the length of the hamming window used to compute the spectrogram 

params.SpectrogramOverlap = 1;                   % number of points of overlap for the spectrogram computation (entire period)

params.EventSpectrogramOverlap = 0;              % number of points of overlap for the spectrogram computation (event period)

params.StimulationArtifactPeriod = 10;           % time in seconds of the period after stimulation within which to remove stimulation artifacts 

params.SeizureOnsetTimeSkip = 2.5;               % time in seconds to skip after seizure onset (bc of poor signal quality right after seizure onset)

params.SpectrogramFrequencyRange = [0 100];      % frequency (hz) range over which to calculate the time course spectrogram 

params.nfft = 2048;                              % number of fft points to evaluate the fft for event analysis 

params.SpectrogramNfft = 8000;                   % number of fft points to evalue the spectrogram function for event analysis 

params.ClickInterval = 4;                        % minimum time interval in second between two clicks, the earlier click will be discarded if two clicks are too close 

params.SpectLogscale = true;                     % liner (watt/Hz) or log scale (db/Hz) for the spectrograms, if true then 10*log10() will be taken on the PSD     

%% Event response processing 

params.EventResponseProcessing.nDeviationsAboveMean = Inf;     % how many stds above mean to use for click exclusion

params.EventResponseProcessing.percentOutlierValues = 1;       % percentage of outlier values to use for click exclusion 

params.EventResponseProcessing.freqRangeOfInterest = [0, 2];   % frequency range over which outlier values will be used to exclude clicks 

%% Normalization
params.normalize = false;                   % whether or not to normalize the signals for spectrograms (TC and click-based)

params.Normalization.Method = 'seizure';        % the method used to normalize signals by baseline. possible values are ['seizure', 'animal']

params.Normalization.ClickBasedMethod = 'seizure';        % the method used to normalize signals by baseline. possible values are ['seizure', 'animal']

params.Normalization.Zero = true;           % whether or not to normalize offset as well

%% Analysis 
params.Analysis.ByAnimalFirst = true;           % average over animals individually and then perform the analysis on those averages 

params.Analysis.ByAnimalValue = 'average';      % whether to average first or take the median first. possible values are ['average', 'median']

params.Analysis.summarystats = 'average';       % EP summary statistics ['average', 'median']

params.Analysis.errorbar = 'std';               % EP summary statistics ['std', 'sem']

params.Analysis.BehPlotStats = 'average';       % Behavior plots whether to show boxplots (median) or bars (average) ['average', 'median']

params.Analysis.BehPlotError = 'sem';           % Behavior plots ['std', 'sem']

params.Analysis.MUAstats = 'average';           % MUA summary statistics ['average', 'median']

%% save 
save(fullfile(savePath, 'params.mat'), 'params'); 
end