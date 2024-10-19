% Get paramters
% Edited by Jerry Liu, 09/20/2022
function Parameter = GetParameter_220919(HCstim,savePath)
    % this function get the suitable parameters for analysis. call this
    % function before implementing the analysis and pass the output struct 
    % as an argument to functions that need it. if you are to change any of
    % parameters, change them here, not within functions. The only optional
    % argument here is the save path to which to save the params mat file. 
    % If a save path is not specified, will save to the MATLAB path. 
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
    %   > IctalCutoff
    %   > ContrlOffset
    % Created by Abdo Sharaf - abdo.sharaf@yale.edu, July 25th 2020
    % Xinyuan Zheng - edit June 2 2022
    % *********************************************************************
    %% Initialization
    % Sampling rate.sps
    Parameter.SampleRate = 1000; % LFP:1k, MUA:20k
    % Regard clicks in baseline as "hit" or not (or reference clicks really)
    Parameter.BaselineRef = true; % true: all 'a' clicks are defined as "hit"                  
    % Duration of the pre-ictal (baseline) period 
    Parameter.PreIctalOffset = 60; %[s]
    % Duration of the post-ictal period
    Parameter.PostIctalOffset = 30; %[s]
    % Ictal Cutoff
    Parameter.IctalCutoff = 15; %[s]
    % Duration of inter-mediate state (delay) between post-ictal and recovery 
    Parameter.RecoveryDelay = 0; %[s]                       
    % Duration of the recovery period 
    Parameter.RecoveryPeriod = 60; %[s] 
    % Duration of control period, same with ictal cutoff as reference
    Parameter.ControlOffset = 15; %[s]
    % Checking ipsi and contral, choose an end-time for Behavior-Analysis
    % 'ipsi': szEnd3
    % 'contra': szEnd4 
    % 'max_ipsi_contra': max([szEnd1,szEnd2,szEnd3,szEnd4])
    if strcmp(HCstim,'ipsi')
        % Electrode locations of interest, other locations are discarded
        % 'ipsiLHC', 'ipsiRHC', 'ipsi_all', 'contra', 'all'
        Parameter.electrodeLocation = 'ipsi_all';
        Parameter.endTimeMethod = 'ipsi'; 
    elseif strcmp(HCstim,'contra')  
        Parameter.electrodeLocation = 'contra';
        Parameter.endTimeMethod = 'contra';
    elseif strcmp(HCstim,'all')
        Parameter.electrodeLocation = 'all';
        Parameter.endTimeMethod = 'max_ipsi_contra';
    else
        disp('Wrong Stimulation Input');
    end
    % First lick delay threshold for classifying using "first_lick" method
    Parameter.firstLickWindow = 1.05; %[s]                   
    % Time window before and after "a" click for calculating max lick rate 
    Parameter.maxLickRateWindow = [-0.1, 1.05]; %[s]        
    % Bin size used for estimating lick rate with binning 
    Parameter.binSize = 0.1; %[s]
    % Offset before "a" click for Behavior-Analysis
    Parameter.PreStimOffset = 2; %[s]                     
    % Offset after "a" click for Behavior-Analysis
    Parameter.PostStimOffset = 2+0.1; %[s]                      
    % Offset before "a" click for Electrophysiology-Analysis
    Parameter.EPPreStimOffset = 2; %[s]                     
    % Offset after "a" click for Electrophysiology-Analysis
    Parameter.EPPostStimOffset = 2; %[s]
    % Time window for "a" click beyond which responses are "timed out" 
    Parameter.trialWindow = 3; %[s]                             
    % Temporal resolution in the process of estimating lick rate 
    Parameter.TimeCourseResolution = 0.01; %[s]             
    % Sliding window size estimating lick rate in 'sliding window' method 
    Parameter.SlidingWindowSize = 0.1; %[s]
    % Lick rate threshold for classifiying in 'window lick rate' method
    Parameter.WinLickRateThreshold = 0.5;                
    % Click classification method:
    % 'first lick'; 'window lick rat'; 'sliding window' 
    Parameter.ClassificationMethod = 'first_lick';      
    % Baseline fraction when classifying with normalized lick rate
    Parameter.FractionOfBaseline = 0.5;                 
    % Hamming window lenght for spectrogram 
    Parameter.SpectrogramWindow = 1000;                 
    % Number of overlap points for spectrogram (entire period)
    Parameter.SpectrogramOverlap = 1;                   
    % Number of overlap points for spectrogram (event period)
    Parameter.EventSpectrogramOverlap = 0; 
    % Time after stimulation within which to polyfit stimulation artifacts 
    % (poor signal quality after sz onset)
    Parameter.StimulationArtifactPeriod = 10; %[s]
    % Time skip after sz induction onset
    Parameter.SeizureOnsetTimeSkip = 2.5;  %[s]             
    % Frequency range for the time course spectrogram
    Parameter.SpectrogramFrequencyRange = [1 50]; % [Hz]      
    % Number of fft points for event analysis 
    Parameter.nfft = 2048;                              
    % Number of fft points for spectrogram function in event analysis 
    Parameter.SpectrogramNfft = 8000;                   
    %% Event response processing 
    % stds threshold above mean to use for click exclusion
    Parameter.EventResponseProcessing.nDeviationsAboveMean = Inf;     
    % Percentage of outlier values to use for click exclusion 
    Parameter.EventResponseProcessing.percentOutlierValues = 1;
    % Frequency range over which outlier values used to exclude "a" clicks
    Parameter.EventResponseProcessing.freqRangeOfInterest = [0,2];   
    %% Normalization
    % Normalization Method to normalize signals by baseline. 
    Parameter.Normalization.Method = 'seizure';  % 'seizure' or 'animal'
    % Whether or not to normalize offset
    Parameter.Normalization.Zero = false;   % true: normalize offset
    %% Analysis 
    % average over animals individually and then analyze averages
    Parameter.Analysis.ByAnimalFirst = false;    
    % Average first or take the median first
    Parameter.Analysis.ByAnimalValue = 'average'; % 'average', 'median'
    %% .mat save 
    save(fullfile(savePath,'Parameter.mat'),'Parameter'); 
end