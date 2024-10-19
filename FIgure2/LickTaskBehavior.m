function Behavior = LickTaskBehavior(seizurePaths, seizureCell, params, varargin)
% this function extracts all the behavioral infromation used later for
% analysis. The notation [n] refers to the length n of arrays when
% applicable.
% INPUTS:
%   - seizurePaths [nSeizures]: a string array with the full paths to the
%   seizure files. This is one of the outputs of the getAllSeizureInfo
%   function.
%   - seizureCell: one of the outputs of getAllSeizureInfo function. See
%   its documentation for more info.
%   - params: a struct that contains *at least* the following fields:
%       > firstLickWindow: a time limit (in seconds) on the first lick within
%        which a click is classified as a hit. If the first lick is beyond
%        this window, the click event is classified as a miss
%       > maxLickRateWindow: this is an array of the form [st, nd] where st
%       is the time relative to the onset of the click that indicates the
%       start of the window that will be used to calculate max lick rate
%       associated with a click, and nd indicates the end of that window
%       relative to the click onset. This will be used for click
%       classification using the maximum lick rate method and more
%       generally to calculate the maximum lick rate for a certain period
%       in a certain seizure file.
%       > binSize: the size of the bin/interval (in seconds) which will be
%       used to calculate lick rate timecourses
%       > preStimOffset: time (in seconds) before the click event to
%       include in the analysis window
%       > postStimOffset: time (in seconds) after the click event to
%       include in the analysis window
%       > Note: this struct is what will be passed onto the
%       getSeizureTimingInfo function. See the documentation for that
%       function for additional fields that this strucut should contain.
%   - savePath: this is an optional argument that will specify where to
%   save the Behavior structure. The structure will be save under the name
%   'Behavior.mat'. If not specified, the file will be save in the MATLAB
%   path.
% OUTPUTS:
%   -   Behavior: a struct with all relevant behavioral information used
%   for analysis. This struct is utilitzed by the AnalyzeBehavior function.
% Created By: Abdo Sharaf - abdo.sharaf@yale.edu
% Modified By: Xinyuan Zheng - xinyuan.zheng@yale.edu
%**************************************************************************

%% Parse in optional arguments
savePath = pwd;
if nargin > 3
    savePath = varargin{1};
end

%% Define and initialize variables

% parse parameters from the params input struct
maxLickRateWindow = params.maxLickRateWindow;
binSize = params.binSize;
preStimOffset = params.PreStimOffset;
postStimOffset = params.PostStimOffset;
clickInterval = params.ClickInterval;

% include/exclude array
toInclude = string(seizureCell(2:end,6));

% array of time points indicating the start time (sec) of each bin
% binTimePts = -preStimOffset : binSize : postStimOffset;
binTimePts = -preStimOffset : binSize : postStimOffset-binSize;

% number of bins (i.e. width of the analysis window in bins)
nBins = length(binTimePts);

% cell array for all seizure files that holds delay to first lick for each
% period (each an entry is an array with length = no. of hit/miss clicks
% for that period)
FirstLickAll.Miss = cell(length(seizurePaths), 4); % delay to first lick for the misses
FirstLickAll.Hit = cell(length(seizurePaths), 4);  % delay to first lick for the hits
FirstLickAll.AllClicks = [];    % delay to the first lick for all clicks

% same as above for the total number of licks
NumLicksAll.Miss = cell(length(seizurePaths), 4); % total number of licks for the misses
NumLicksAll.Hit = cell(length(seizurePaths), 4);  % total number of licks for the hits
NumLicksAll.AllClicks = [];     % total number of licks for all the clicks

% this is a (nClicks x 3) array where nClicks is the total number of seizure
% clicks in all seizure files. Each row corresponding to a click.
% the columns correspond to: 1. the index of the seizure file from which
% click comes, 2. the time from seizure, 3. whether the click is a hit or a miss
TimeFromSzStartAll = [];

% array to hold the duration of all the included seizures
IncludedSzDuration = [];

% array to hold the click times for all seizures and periods
allClickTimes = [];

% time courses for hits and misses. these are time courses of the number of licks
% this is 3D array where the first dimension refers to the bin number,
% the second refers to the seizure file and the third refers to the period
TimeCourse.Hit = nan(nBins, sum(strcmp(toInclude,'Include')), 4);   % avg across hits only
TimeCourse.Miss = nan(nBins, sum(strcmp(toInclude,'Include')), 4);  % avg across misses only
TimeCourse.All = nan(nBins, sum(strcmp(toInclude,'Include')), 4); % this will be average across hits and misses

% time courses for individual clicks.
ClickTimeCourses.Hit = [];
ClickTimeCourses.Miss = [];
ClickTimeCourses.EstimatedLickRate = {};
ClickTimeCourses.AllClicks = [];

% array to store the maximum lick rate for each period for each seizure
% file. There are two ways to look at lick rate: (1) look at the maximum
% rate within a window (specified in parameters), this rate is basically
% the maximum number of licks divided by bin size, (2) look at the overall
% lick rate between the time of the first lick and the time of the next
% click, this will be the total number of licks divided by that time
% difference.
LickRate.Max.PerSeizure = nan(sum(strcmp(toInclude,'Include')), 4);
LickRate.Max.PerClick = [];
LickRate.Overall.PerSeizure = nan(sum(strcmp(toInclude,'Include')), 4);
LickRate.Overall.PerClick = [];
LickRate.NormWindowAvg = [];

% wheel speed
WheelSpeed.Max.PerClick = [];
WheelSpeed.Mean.PerClick = [];
WheelSpeed.AllClicks = {};

% array to store the average of intervals between clicks in a certain
% period for each seizure files
avgClickIntervals = nan(sum(strcmp(toInclude,'Include')), 4);

% array to store the total number of clicks in each period for each seizure
% file
NumClicksAll = zeros(sum(strcmp(toInclude,'Include')), 4);

% array to store the indices of the included seizure files
seizures = [];

% array to store the names of animals to which each seizure belongs
Animal = string(seizureCell(2:end,1))';

% this array will hold a summary of click information for each period for
% all seizures. The number of columns is (nPeriods * 5 + 1) and the number
% of rows is the total number of included seizures. The first column hold
% the index of the seizure file. Each five columns after the first column
% will then hold info about each period (i.e. the first 5
% cols after column 1 will be for baseline, the second 5 cols will be for
% ictal etc.). The info held in the five cols are (1) avg delay to first
% lick for the hits, (2) number of hits, (3) avg delay to first lick for the
% misses, (4) number of misses, (5) percent hits (or hit rate)
ClickResponses = zeros(sum(strcmp(toInclude,'Include')), 26);

    
%% Behavioral segmentation and pre-processing
% for each seizure file, get data associated with clicks from each period
% (i.e baseline, ictal, postictal, recovery, and control)

szCount = 0; % counter for the number of included seizures
for file = 1:length(seizurePaths)
    % file = 69

    % path to the seizure file
    % filepath = 'F:\xinyuan\Mouse Project\Behavior_partial_confirmed_LOVO\Tomoko\200723_MUA_Tomoko_000.mat'
    filepath = seizurePaths{file};
    disp(filepath);
    % check if this file should be included or not. if not, continue to
    % next iteration
    if strcmp(toInclude(file), 'Include')
        szCount = szCount + 1;  % increment the number of seizures
    else
        continue;   % move on to next iteration
    end
    
    % get the click, lick, and timing information for this seizure file
    seizureInfo = getSeizureTimingInfo(filepath, seizureCell, params);
    
    % parse useful variables from the seizureInfo structure
    AllTimes = seizureInfo.AllTimes;
    Lick = seizureInfo.Lick;
    Click = seizureInfo.Click;
    ClickInMs = seizureInfo.ClickInMs;
    if seizureInfo.hasWheel
        wheelData = seizureInfo.Wheel;
    end
    
    % store the index for this seizure
    seizures = [seizures, file];
    
    % update the duration of this seizure
    IncludedSzDuration = [IncludedSzDuration; [file,...
        ((AllTimes(4)-AllTimes(3))/params.SampleRate)]];
    
    % loop through the periods
    for period = 1:5 % period = 4
        % period = 5
        % check if there is a control period. if not continue
        if ~seizureInfo.hasControl && period == 5
            continue;
        end
        
        % logical index for the clicks within this period's time range
        IncludedClicks = ClickInMs > AllTimes(1, 2*period-1) & ...
            ClickInMs < AllTimes(1, 2*period);
        
        % get click times for the included clicks
        % [for now, we are not checking for duplicate click events. this is
        % something that we may wanna implement in the future?]
        clickArray = Click.times(IncludedClicks);
        
        % [fixed the above question. xinyuan 25 Jul 2022]
        % remove dulicated clicks
        clickArray(diff(clickArray) < clickInterval) = []; 
        
        
        % number of clicks and click intervals
        if isempty(clickArray) % if no clicks in this period
            numClicks = 0;
            afterLastClick = Lick.times(end);
        else
            % get the time of the click after the last click (for calculating
            % number of clicks and click intervals
            if find(Click.times==clickArray(end))+1 <= length(Click.times)
                afterLastClick = Click.times(find(Click.times==clickArray(end))+1);
            else
                afterLastClick = Lick.times(end);   % if it's the last click then the upper bound will be the end of the trial (last lick time)
            end
            numClicks = numel(clickArray);
            clickIntervals = diff([clickArray; afterLastClick]);
            avgClickIntervals(szCount, period) = nanmean(clickIntervals);
        end
        
        if numClicks ~= 0
            % classify and analyze click events
            if ~isempty(Lick.times) % only if the seizure file actually has lick times
                switch params.ClassificationMethod
                    case 'first_lick'
                        if ~params.BaselineRef
                            %analyze all clicks and classify based on first_licks
                            [allClickClasses, responseTimes, allClickNumLicks, allClickTimeCourses] = analyzeClicks(clickArray,...
                                afterLastClick, Lick.times, binTimePts, params, params.ClassificationMethod);
                        else
                            if period == 1
                                allClickClasses = ones(length(clickArray), 1);  % all clicks in baseline are regarded as hits
                                % we are just getting the other info. we don't really care about
                                % classification with baseline in this case
                                [~, responseTimes, allClickNumLicks, allClickTimeCourses] = analyzeClicks(clickArray,...
                                    afterLastClick, Lick.times, binTimePts, params, 'first_lick');
                            else
                                %analyze all clicks and classify based on first_licks
                                [allClickClasses, responseTimes, allClickNumLicks, allClickTimeCourses] = analyzeClicks(clickArray,...
                                    afterLastClick, Lick.times, binTimePts, params, params.ClassificationMethod);
                                
%                                 if period == 2
%                                         SzBeh.seizureID{file} = seizureCell{file+1,2};
%                                         SzBeh.ClickType{file} = allClickClasses;
%                                         SzBeh.Impaired{file} = sum(allClickClasses == 0 | allClickClasses == -1)/length(allClickClasses);
%                                         SzBeh.Spared{file} = sum(allClickClasses == 1)/length(allClickClasses);
%                                 end
                            end
                        end
                        
                    case 'win_lick_rate_and_first_lick'
                        % analyze all clicks and classify based on both window lick
                        % rate and delay to the first lick
                        [allClickClasses, features, allClickNumLicks, allClickTimeCourses] = analyzeClicks(clickArray,...
                            afterLastClick, Lick.times, binTimePts, params, params.ClassificationMethod);
                        
                        responseTimes = features(:, 1);
                        
                    case 'norm_lick_rate_and_first_lick'
                        % analyze all clicks and classify based on both the
                        % normalized window lick rate and delay to first lick
                        % (BaselineRef has to be true here since we are normalizing by baseline)
                        if period == 1  % baseline gets special treatment here
                            allClickClasses = ones(length(clickArray), 1); % all clicks in baseline are regarded as hits
                            
                            % we'll use the analyzeClicks function with a
                            % different classification method just cuz we want
                            % the other info. we don't really care about
                            % classification with baseline in this case
                            [~, responseTimes, allClickNumLicks, allClickTimeCourses] = analyzeClicks(clickArray,...
                                afterLastClick, Lick.times, binTimePts, params, 'first_lick');
                            
                            % define the timecourses array that will be used
                            % for normalization
                            normBy = allClickTimeCourses;
                            
                            % get the norm win lick rate
                            stBinIdx = ceil(preStimOffset / binSize);
                            pdWinLickRate = sum(allClickTimeCourses(:, stBinIdx:end), 2)/postStimOffset;
                            bslAvg = nanmean(pdWinLickRate);
                            normWinRate = pdWinLickRate./bslAvg;
                            
                        else
                            [allClickClasses, features, allClickNumLicks, allClickTimeCourses] = analyzeClicks(clickArray,...
                                afterLastClick, Lick.times, binTimePts, params, params.ClassificationMethod, normBy);
                            
                            responseTimes = features(:, 1);
                            
                            % get the norm win lick rate for
                            normWinRate = features(:, 2);
                        end
                        
                    otherwise
                        error('Invalid Classification Method Specified');
                end
                
                % estimate the lick rate time course with a sliding window
                ClickTimeCourses.EstimatedLickRate = [ClickTimeCourses.EstimatedLickRate;...
                    estimateLickRate(clickArray, afterLastClick, Lick.times, params, 'sliding_window')];
                
                
                % get indices for misses and hits
                missIndcs = allClickClasses == 0 | allClickClasses == -1;
                hitIndcs = allClickClasses == 1;
                
                %%% initialize the arrays that will hold this iteration's information
                % time to first lick
                hitFirstLick = responseTimes(hitIndcs);  % collects time to first lick for hit trials in this period
                missFirstLick = responseTimes(missIndcs); % collects time to first lick for miss trials in thei period
                
                % number of licks
                hitNumLicks = allClickNumLicks(hitIndcs);   % collects the number of licks for hit trials in this period
                missNumLicks = allClickNumLicks(missIndcs);  % collects the number of licks for miss trials in this period
                
                %%% all timecourses
                % collects all timecourses for hits in this period (i.e. nrows = nhits)
                hitTimeCourses = allClickTimeCourses(hitIndcs, :);
                % collects all timecourses for misses in this period
                missTimeCourses = allClickTimeCourses(missIndcs, :);
                
                % append to the clicktimecourse array
                ClickTimeCourses.Hit = [ClickTimeCourses.Hit;...
                    [hitTimeCourses, file*ones(size(hitTimeCourses,1), 1),...
                    period*ones(size(hitTimeCourses,1), 1)]];
                ClickTimeCourses.Miss = [ClickTimeCourses.Miss;...
                    [missTimeCourses, file*ones(size(missTimeCourses,1), 1),...
                    period*ones(size(missTimeCourses,1), 1)]];
                ClickTimeCourses.AllClicks = [ClickTimeCourses.AllClicks;...
                    [allClickTimeCourses, file*ones(numClicks, 1),...
                    period*ones(numClicks, 1), allClickClasses]];
                
                % modify the classes array to replace -1 with 0
                MissOrHit = allClickClasses;
                MissOrHit(MissOrHit==-1) = 0;
                
                % calculate and store the max lick rate per click
                allClickMaxLickRate = max(allClickTimeCourses(:, binTimePts >= maxLickRateWindow(1)...
                    & binTimePts <= maxLickRateWindow(2)),[], 2)./binSize;
                thisMaxLickRateByClick = [file*ones(numClicks,1),...
                    allClickMaxLickRate, period*ones(numClicks,1), MissOrHit];
                LickRate.Max.PerClick = [LickRate.Max.PerClick; thisMaxLickRateByClick];
                
                % calculate and store the overall lick rate for each click
                allClickOverallLickRate = allClickNumLicks ./ (clickIntervals - responseTimes);
                allClickOverallLickRate(isnan(allClickOverallLickRate)) = 0;   % nan lick rate are actually zero since they have no associated licks
                LickRate.Overall.PerClick = [LickRate.Overall.PerClick;...
                    [file*ones(numClicks,1), allClickOverallLickRate, period*ones(numClicks,1),...
                    MissOrHit]];
                
                % store the normalized window avg lick rate
                if strcmp(params.ClassificationMethod, 'norm_lick_rate_and_first_lick')
                    LickRate.NormWindowAvg = [LickRate.NormWindowAvg;...
                        [file*ones(numClicks,1), normWinRate, period*ones(numClicks,1),...
                        MissOrHit]];
                end
                
                % store the click times
                allClickTimes = [allClickTimes; [file*ones(numClicks,1),...
                    clickArray, period*ones(numClicks,1), MissOrHit]];
                
                % wheel speed
                if seizureInfo.hasWheel
                    clickWheelData = sequenceWheelData(clickArray, wheelData, params);
                    % smooth out aburpt changes
                    clickWheelData_smoothed = cellfun(@(v)smoothdata(v, 'sgolay'),...
                        clickWheelData, 'UniformOutput', false);
                    % calculate the speed
                    clickWheelSpeed = cellfun(@(v)abs(movingslope(v, 5, 2)), clickWheelData_smoothed, 'UniformOutput', false);
                    % remove outliers (maybe not necessary if we smooth?)
                    clickWheelSpeed = cellfun(@(v)filloutliers(v, 'nearest'), clickWheelSpeed, 'UniformOutput', false);
                    %                 for clk = 1:numClicks
                    %                     % get peaks
                    %                     [~, locs] = findpeaks(clickWheelSpeed{clk},'MinPeakProminence', 5);
                    %                     % remove peaks
                    %                     clickWheelSpeed{clk}(locs) = NaN;
                    %                     % fill missing points
                    %                     clickWheelSpeed{clk} = fillmissing(clickWheelSpeed{clk}, 'nearest');
                    %                 end
                    WheelSpeed.AllClicks = [WheelSpeed.AllClicks; clickWheelSpeed];
                    clickWheelSpeedAvg = cellfun(@nanmean, clickWheelSpeed);
                    clickWheelSpeedMax = cellfun(@nanmax, clickWheelSpeed);
                    
                else
                    WheelSpeed.AllClicks = [WheelSpeed.AllClicks; num2cell(nan(numClicks, 1))];
                    [clickWheelSpeedAvg, clickWheelSpeedMax] = deal(nan(numClicks, 1));
                end
                
                WheelSpeed.Max.PerClick = [WheelSpeed.Max.PerClick;...
                    [file*ones(numClicks,1), clickWheelSpeedMax, period*ones(numClicks,1),...
                    MissOrHit]];
                WheelSpeed.Mean.PerClick = [WheelSpeed.Mean.PerClick;...
                    [file*ones(numClicks,1), clickWheelSpeedAvg, period*ones(numClicks,1),...
                    MissOrHit]];
                
                % array to store the time from seizure onset for the seizure priod
                % the 1st column is time from seizure onset and the 2nd is whether
                % a click is a hit(1) or a miss(0)
                if period == 2  % if this is a seizure period
                    
                    % calculate the time from seizure onset for all clicks and
                    % append to the All array
                    timeFromSzStart = [file*ones(numClicks,1),...
                        (clickArray*params.SampleRate - (AllTimes(3)))./params.SampleRate,...
                        MissOrHit];
                    TimeFromSzStartAll = [TimeFromSzStartAll; timeFromSzStart];
                end
                
                % number of hits and misses
                hitCount = sum(hitIndcs);
                missCount = sum(missIndcs);
            end
            
            % append to the global arrays (the ones that hold data for all
            % files)
            FirstLickAll.Miss{file, period} = missFirstLick;
            FirstLickAll.Hit{file, period} = hitFirstLick;
            FirstLickAll.AllClicks = [FirstLickAll.AllClicks;...
                [file*ones(numClicks,1), responseTimes, period*ones(numClicks,1), MissOrHit]];
            
            NumLicksAll.Miss{file, period} = missNumLicks;
            NumLicksAll.Hit{file, period} = hitNumLicks;
            NumLicksAll.AllClicks = [NumLicksAll.AllClicks;...
                [file*ones(numClicks,1), allClickNumLicks, period*ones(numClicks,1), MissOrHit]];
            
            ClickResponses(szCount, 1) = file;
            ClickResponses(szCount, (2:6)+(5*(period-1))) = [nanmean(hitFirstLick),...
                hitCount, nanmean(missFirstLick), missCount, (hitCount/(hitCount + missCount))];
            
            TimeCourse.Hit(:, szCount, period) = nanmean(hitTimeCourses,1)';
            TimeCourse.Miss(:, szCount, period) = nanmean(missTimeCourses,1)';
            TimeCourse.All(:, szCount, period) = nanmean([hitTimeCourses; missTimeCourses]);
            
            LickRate.Max.PerSeizure(szCount, period) = max(TimeCourse.All(binTimePts >= maxLickRateWindow(1)...
                & binTimePts <= maxLickRateWindow(2), szCount, period)) / binSize;
            
            LickRate.Overall.PerSeizure(szCount, period) = mean(LickRate.Overall.PerClick((LickRate.Overall.PerClick(:,1)...
                == szCount) & (LickRate.Overall.PerClick(:,3)== period), 2));
            
            NumClicksAll(szCount, period) = numClicks;
        end
    end
    
end

%% Build the Behavior structure returned by the function

%%% Averages by click
Behavior.ByClick.TimeFromSzStart = TimeFromSzStartAll;
Behavior.ByClick.FirstLickAll = FirstLickAll;
Behavior.ByClick.NumLicksAll = NumLicksAll;
Behavior.ByClick.LickRateAll.Max = LickRate.Max.PerClick;
Behavior.ByClick.LickRateAll.Overall = LickRate.Overall.PerClick;
Behavior.ByClick.LickRateAll.NormWindowAvg = LickRate.NormWindowAvg;
Behavior.ByClick.ClickTimes = allClickTimes;
Behavior.ByClick.TimeCourse = ClickTimeCourses;
if ~isempty(allClickTimes)
    Behavior.ByClick.AnimalsAll = Animal(allClickTimes(:, 1));
else
    Behavior.ByClick.AnimalsAll = [];
end
Behavior.ByClick.WheelSpeed.Max = WheelSpeed.Max.PerClick;
Behavior.ByClick.WheelSpeed.Mean = WheelSpeed.Mean.PerClick;
Behavior.ByClick.WheelSpeed.All = WheelSpeed.AllClicks;

%%% Averages by seizures

% note that this array will now have 4 cols where the first col corresponds
% to seizure file index, the second col to the number of clicks in the
% seizure period of that file, the third col to the average time from
% seizure onset across clicks for that seizure, and the 4th col to the hit
% rate (basically percentage of hit clicks) for the seizuer period in that
% file
Behavior.BySeizure.TimeFromSzStart = getAvgByClass(TimeFromSzStartAll(:,2:end), TimeFromSzStartAll(:,1));

Behavior.BySeizure.IncludedSzDuration = IncludedSzDuration;
Behavior.BySeizure.TimeCourse = TimeCourse;

Behavior.BySeizure.MeanTimeCourse.Hit = [nanmean(TimeCourse.Hit(:, :, 1), 2),...
    nanmean(TimeCourse.Hit(:, :, 2), 2), nanmean(TimeCourse.Hit(:, :, 3), 2),...
    nanmean(TimeCourse.Hit(:, :, 4), 2)];
Behavior.BySeizure.MeanTimeCourse.Miss = [nanmean(TimeCourse.Miss(:, :, 1), 2),...
    nanmean(TimeCourse.Miss(:, :, 2), 2), nanmean(TimeCourse.Miss(:, :, 3), 2),...
    nanmean(TimeCourse.Miss(:, :, 4), 2)];

Behavior.BySeizure.HitRate = ClickResponses(:, 6:5:21);
Behavior.BySeizure.HitLatency = ClickResponses(:,2:5:17);
Behavior.BySeizure.MissLatency = ClickResponses(:,4:5:19);
Behavior.BySeizure.DurationAndHitRate = [ClickResponses(:, 1), IncludedSzDuration, ClickResponses(:, 11)];

Behavior.BySeizure.ClickIntervals = avgClickIntervals;
Behavior.BySeizure.NumClicks = NumClicksAll;

Behavior.BySeizure.LickRateAll.Max = LickRate.Max.PerSeizure;
Behavior.BySeizure.LickRateAll.Overall = LickRate.Overall.PerSeizure;


%%% Averages by animal
Behavior.ByAnimal.IncludedSzDuration = getAvgByClass(IncludedSzDuration(:,2), Animal(IncludedSzDuration(:, 1)));
Behavior.ByAnimal.HitRate = getAvgByClass(ClickResponses(:,6:5:21), Animal(seizures'));
Behavior.ByAnimal.HitLatency = getAvgByClass(ClickResponses(:,2:5:17), Animal(seizures'));
Behavior.ByAnimal.DurationAndHitRate = getAvgByClass([IncludedSzDuration(:, 2),...
    ClickResponses(:, 11)], Animal(ClickResponses(:,1)));

% note that these time courses now will be of dimensions (nBins+2 x
% nAnimals) where the first two row refer to animal and total number of clicks
% associated with that animal for a certain period, respectively. Row
% 2:nBins+2 will correspond to the actual lick rate time course. Each
% column will correspond to one animal
Behavior.ByAnimal.TimeCourse.Hit.Baseline = getAvgByClass(TimeCourse.Hit(:, :, 1)', Animal(seizures)')';
Behavior.ByAnimal.TimeCourse.Hit.Ictal = getAvgByClass(TimeCourse.Hit(:, :, 2)', Animal(seizures)')';
Behavior.ByAnimal.TimeCourse.Hit.PostIctal = getAvgByClass(TimeCourse.Hit(:, :, 3)', Animal(seizures)')';
Behavior.ByAnimal.TimeCourse.Hit.Recovery = getAvgByClass(TimeCourse.Hit(:, :, 4)', Animal(seizures)')';

Behavior.ByAnimal.TimeCourse.Miss.Baseline = getAvgByClass(TimeCourse.Miss(:, :, 1)', Animal(seizures)')';
Behavior.ByAnimal.TimeCourse.Miss.Ictal = getAvgByClass(TimeCourse.Miss(:, :, 2)', Animal(seizures)')';
Behavior.ByAnimal.TimeCourse.Miss.PostIctal = getAvgByClass(TimeCourse.Miss(:, :, 3)', Animal(seizures)')';
Behavior.ByAnimal.TimeCourse.Miss.Recovery = getAvgByClass(TimeCourse.Miss(:, :, 4)', Animal(seizures)')';

Behavior.ByAnimal.NumClicks = getAvgByClass(NumClicksAll, Animal(seizures));

% save
save(fullfile(savePath, 'Behavior.mat'),'Behavior');

end

%% Helper functions
function avgs = getAvgByClass(A, classIDs)
% this function averages an array by 'type'. For example, class could be
% seizure, in which case the array elements will correspond to certain
% seizures and all elements corresponding to one seizure will be averaged.
% Same thing if the class was animal.
% Inputs:
%   - A: array to be averaged
%   - classIDs: identifiers for each element in A which will be used to
%   group elements in the same class together
% Outputs:
%   - avgs: averaged array across elements of the same classes
%**************************************************************************

% check if the inputs are the same size
if size(A, 1) ~= length(classIDs)
    error('A and classIDs must be the same size');
end

% initialize ouputs
avgs = [];

% get the uniqueIDs for the classes
% if the classIDs is a string array, get numeric IDs for each class.
% A class's ID will be the index of the first occurence of that class in the classIDs array
if isstring(classIDs)
    [~, unqIDs, ~] = unique(classIDs);
else
    unqIDs = unique(classIDs);
end

for i = 1:length(unqIDs)
    avgs(i, 1) = unqIDs(i);     % ID corresponding to this average
    % number of elements in this class
    if isstring(classIDs)
        avgs(i, 2) = sum(classIDs == classIDs(unqIDs(i)));
        avgs(i, 3:2+size(A,2)) = nanmean(A(classIDs == classIDs(unqIDs(i)), :)); % average
        
    else
        avgs(i, 2) = sum(classIDs == unqIDs(i));
        avgs(i, 3:2+size(A,2)) = nanmean(A(classIDs == unqIDs(i), :)); % average
    end
end

end
