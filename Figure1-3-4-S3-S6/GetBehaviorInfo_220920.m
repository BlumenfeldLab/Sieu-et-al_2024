% Get behavior related information
% Edit by Jerry. 09/20/2022
function Behavior = GetBehaviorInfo_220920...
    (szCell,szDir,szTiming,params,SaveDir)
    %*********************************************************************
    % this function extracts all the behavioral information 
    % The notation [n] refers to the length n of arrays when applicable.
    % INPUTS:
    % (1) szDir: [string array] with full paths to the seizure files.
    %            One of the outputs of the GetSeizureInfo_220701.
    % (2) szCell: [cell] output of GetSeizureInfo_220701. 
    % (3) szTiming: [structure] changed from cell output of
    %               GetSeizureInfo_220701.
    % (4) params: struct that contains "at least" the following fields:
    %     > firstLickWindow: a time limit[s] on the first lick within
    %                        which a click is classified as a hit. 
    %                        If the first lick is beyond this window, 
    %                        the click event is classified as a miss.
    %     > maxLickRateWindow: an array of the form [st,nd].  
    %       st: the time relative to the onset of the click that indicates 
    %           the start of the window that will be used to calculate 
    %           max lick rate associated with a click. 
    %       nd: the end of that window relative to the click onset.
    %           This will be used for click classification using the
    %           maximum lick rate method and more generally to calculate
    %           the maximum lick rate for a certain period
    %           in a certain seizure file.
    %     > binSize: the size of the bin/interval[s] which will be
    %                used to calculate lick rate timecourses
    %     > preStimOffset: time[s] before the click in the analysis window
    %     > postStimOffset: time[s] after the click in the analysis window
    % (5) SaveDir: specify where to save the Behavior structure
    %*********************************************************************
    % OUTPUTS:
    % (1) Behavior: a struct with all relevant behavioral information used
    %   for analysis. This struct is used by the AnalyzeBehavior function.
    %     > byClick: done
    %     > bySeizure: to do
    %     > byAnimal: do not care
    % Created By: Abdo Sharaf - abdo.sharaf@yale.edu
    %*********************************************************************
    %% Define and initialize variables
    % 4 refers to 4 periods: baseline,ictal,post-ictal and recovery 
    % parse parameters from the input params struct
    maxLickRateWindow = params.maxLickRateWindow; % [-0.1, 1.05] s
    binSize = params.binSize; % 0.1s
    preStimOffset = params.PreStimOffset; % 2s
    postStimOffset = params.PostStimOffset; % 2s+0.1s
    % include/exclude array
    toInclude = string(szCell(2:end,3));
    % Time points array indicating the start time [s] of each bin
    binTimePts = -preStimOffset:binSize:postStimOffset-binSize;
    % Number of bins (i.e. width of the analysis window in bins)
    nBins = length(binTimePts);
    % cell array for first lick delay: [sz#, period#]. (ByClick)
    % for each sz (row) and each period (column).
    % each entry is an array [hit/miss clicks #] for 1 of 4 periods
    FirstLickAll.Miss = cell(length(szDir),4);% [cell] for miss clicks(>1s)
    FirstLickAll.Hit = cell(length(szDir),4); % [cell] for hit clicks(<1s)
    FirstLickAll.AllClicks = [];    % [array] for all clicks
    % the total number of licks w.r.t the click ahead. (ByClick)
    NumLicksAll.Miss = cell(length(szDir),4); % [cell] for miss clicks(>1s)
    NumLicksAll.Hit = cell(length(szDir),4);  % [cell] for hit clicks(<1s)
    NumLicksAll.AllClicks = [];     % [array] for all clicks
    % TimeFromSzStartAll: (nClicksX3) [arrays]. (ByClick)
    % nClicks: 
    % the number of clicks in all sz files, each row is a click.
    % 3 columns: 
    % 1: the seizure file from which click comes.
    % 2: the time delay from seizure onset.
    % 3: the click is a hit or a miss (1/0).
    TimeFromSzStartAll = [];
    % Durations array of all seizures. (ByAnimal)
    IncludedSzDuration = [];
    % Click times for all seizures and periods. (ByClick.ClickTimes)
    allClickTimes = [];
    % Licks number time-course for hits and misses and is a 3D array
    % the first dimension: bins number,
    % the second dimension: sz files number 
    % the third dimension: the period avg across hits only
    % (BySeizure)
    TimeCourse.Hit = nan(nBins,sum(strcmp(toInclude,'Include')),4);   
    TimeCourse.Miss = nan(nBins,sum(strcmp(toInclude,'Include')),4);
    TimeCourse.All = nan(nBins,sum(strcmp(toInclude,'Include')),4); 
    % time courses for individual clicks. (ByClick.TimeCourse)
    ClickTimeCourses.Hit = [];
    ClickTimeCourses.Miss = [];
    ClickTimeCourses.EstimatedLickRate = {}; % sliding window, not used
    ClickTimeCourses.AllClicks = [];
    % The lick rate for each period in each sz file. 
    % There are two types of lick rate: 
    % (1) the lick rate within a trialWindow (specified in parameters): 
    %     the maximum number of licks divided by bin size. (Max)
    % (2) the lick rate between the time of the first lick and next click:
    %     the total licks number divided by the time difference. (Overall)
    % (BySeizure.LickRateAll)
    LickRate.Max.PerSeizure = nan(sum(strcmp(toInclude,'Include')),4);
    LickRate.Overall.PerSeizure = nan(sum(strcmp(toInclude,'Include')),4);
    % (ByClick.LickRateAll)
    LickRate.Max.PerClick = [];    
    LickRate.Overall.PerClick = [];
    LickRate.NormWindowAvg = []; % not used
    % Wheel speed. (ByClick)
    WheelSpeed.Max.PerClick = [];
    WheelSpeed.Mean.PerClick = [];
    WheelSpeed.AllClicks = {};
    % Clicks number in each period in each sz file. (ByClick)
    NumClicksAll = zeros(sum(strcmp(toInclude,'Include')),4);
    % The average of intervals between clicks 
    % in a certain period in each seizure file. (BySeizure.ClickIntervals)
    avgClickIntervals = nan(sum(strcmp(toInclude,'Include')),4);
    % Indices of the seizure files
    seizures = [];
    % Animal names to which each seizure belongs
    Animal = string(szCell(2:end,1))';
    % Summary of click information for each period for all seizures.
    % Column number is (1+nPeriods*5). First column: index of the sz file.
    % Row number is the number of seizures. 
    % After the first column, each five columns hold info about each period
    % (i.e. the first 5 cols after column 1 will be for baseline,
    % the second 5 cols will be for seizure #1, the last 5 for control). 
    % The info in each five cols are: 
    % (1) avg time delay of (first lick to hits).  (2) number of hits.
    % (3) avg time delay of (first lick to misses). (4) number of misses.
    % (5) hits percent (or hit rate)
    ClickResponses = zeros(sum(strcmp(toInclude,'Include')), 26);
    %% Behavioral segmentation and pre-processing
    % For each sz file, get data associated with clicks from each period
    % (i.e baseline, ictal, postictal, recovery, and control)
    szCount = 0; % counter for the number of included seizures
    for fn = 1:length(szDir)
        % Seizure .mat file directory
        filepath = szDir{fn};
        disp(filepath);
        % If this .mat seizure file is included or not. 
        % If not, continue to the next iteration
        if strcmp(toInclude(fn), 'Include')
            szCount = szCount+1;  % increment the number of seizures
        else
            continue;   % move on to next iteration
        end
        % parse useful variables from the szTiming [structure]
        AllTimes = szTiming(szCount,:).AllTimes;
        Lick = szTiming(szCount,:).Lick; % [structure]
        Click = szTiming(szCount,:).Click; % [structure]
        ClickInMs = szTiming(szCount,:).ClickInMs; 
        if szTiming(szCount,:).hasWheel
            wheelData = szTiming(szCount,:).Wheel; % [structure]
        end       
        % store the index for this seizure
        seizures = [seizures,fn]; % [1D array(sz#)]
        % update the duration of this seizure [2D array(sz#,szDur)]
        IncludedSzDuration = [IncludedSzDuration;[fn,...
            ((AllTimes(4)-AllTimes(3))/params.SampleRate)]];
        % loop through periods: baseline,ictal,postictal,recovery,control
        for period = 1:5
            % check if there is a control period (5). if not continue
            if (~szTiming(szCount,:).hasControl) && (period==5)
                continue;
            end
            % logical index for clicks within this period's time range
            IncludedClicks = ClickInMs>AllTimes(1,2*period-1) & ...
                ClickInMs<AllTimes(1,2*period);
            % get all clicks' time for the included clicks
            clickArray = Click.times(IncludedClicks);
            % number of clicks and click intervals
            if isempty(clickArray) % if no clicks in this period
                numClicks = 0;
                afterLastClick = Lick.times(end);
            else
                % get the time of the click after the last click 
                % (for calculating number of clicks and click intervals)
                if (find(Click.times==clickArray(end)) + 1) <= ...
                        (length(Click.times))
                    afterLastClick = Click.times...
                        (find(Click.times==clickArray(end)) + 1);
                else
                    % if it's the last click then the upper bound 
                    % will be the end of the trial (last lick time)
                    afterLastClick=Lick.times(end);   
                end
                numClicks = numel(clickArray);
                % (cat) 'afterLastClick' to the end and then (diff)
                clickIntervals = diff([clickArray;afterLastClick]);
                avgClickIntervals(szCount,period)=...
                    mean(clickIntervals,'omitnan');
            end
            % classify and analyze click events 
            if numClicks ~= 0
                % only if the seizure file actually has lick times
                if ~isempty(Lick.times) 
                    if ~params.BaselineRef
                        [allClickClasses,responseTimes,...
                            allClickNumLicks,allClickTimeCourses]=...
                            AnalyzeClicks_220724(clickArray,...
                            afterLastClick,Lick.times,binTimePts,params);
                    else
                        if period == 1
                            % clicks in baseline are regarded as hits
                            allClickClasses=ones(length(clickArray),1);
                            [~, responseTimes,allClickNumLicks,...
                                allClickTimeCourses] = ...
                                AnalyzeClicks_220724(clickArray,...
                                afterLastClick,Lick.times,...
                                binTimePts,params);
                        else
                            % classify based on first_licks
                            [allClickClasses,responseTimes,...
                                allClickNumLicks,allClickTimeCourses]...
                                = AnalyzeClicks_220724(clickArray,...
                                afterLastClick,Lick.times,...
                                binTimePts,params);
                        end
                    end
                    % estimate the lick rate time course 
                    % with a sliding window (not used) 
                    ClickTimeCourses.EstimatedLickRate = ...
                        [ClickTimeCourses.EstimatedLickRate;...
                        EstimateLickRate_220630(clickArray,...
                        afterLastClick,Lick.times,params)];
                    %
                    % get indices for misses and hits
                    missInds = allClickClasses==0 | allClickClasses==-1;
                    hitInds = allClickClasses==1;
                    %%% initialize arrays for the iteration's information
                    % time to first lick
                    % time to first lick for hit trials in this period
                    hitFirstLick = responseTimes(hitInds);  
                    % time to first lick for miss trials in thei period
                    missFirstLick = responseTimes(missInds); 
                    % number of licks
                    % the number of licks for hit trials in this period
                    hitNumLicks = allClickNumLicks(hitInds);   
                    % the number of licks for miss trials in this period
                    missNumLicks = allClickNumLicks(missInds);
                    %%% all timecourses
                    % hits timecourses in this period (i.e. nrows = nhits)
                    hitTimeCourses = allClickTimeCourses(hitInds,:);
                    % misses timecourses in this period
                    missTimeCourses = allClickTimeCourses(missInds,:);
                    %%% append to the clicktimecourse array
                    ClickTimeCourses.Hit = [ClickTimeCourses.Hit;...
                        [hitTimeCourses,...
                        fn*ones(size(hitTimeCourses,1),1),...
                        period*ones(size(hitTimeCourses,1),1)]];
                    ClickTimeCourses.Miss = [ClickTimeCourses.Miss;...
                        [missTimeCourses,...
                        fn*ones(size(missTimeCourses,1),1),...
                        period*ones(size(missTimeCourses,1),1)]];
                    ClickTimeCourses.AllClicks=...
                        [ClickTimeCourses.AllClicks;...
                        [allClickTimeCourses,fn*ones(numClicks,1),...
                        period*ones(numClicks,1),allClickClasses]];
                    % modify the classes array to replace -1 with 0
                    MissOrHit = allClickClasses;
                    MissOrHit(MissOrHit==-1) = 0;
                    % calculate and store the max lick rate per click
                    allClickMaxLickRate = max(allClickTimeCourses(:,...
                        binTimePts>=maxLickRateWindow(1) & ...
                        binTimePts<=maxLickRateWindow(2)),[],2)./binSize;
                    thisMaxLickRateByClick = [fn*ones(numClicks,1),...
                        allClickMaxLickRate,...
                        period*ones(numClicks,1),MissOrHit];
                    LickRate.Max.PerClick = ...
                        [LickRate.Max.PerClick;thisMaxLickRateByClick];
                    % the overall lick rate for each click
                    allClickOverallLickRate = allClickNumLicks./...
                        (clickIntervals-responseTimes);
                    % nan lick rate are zero (no associated licks)
                    allClickOverallLickRate...
                        (isnan(allClickOverallLickRate))=0;   
                    LickRate.Overall.PerClick=...
                        [LickRate.Overall.PerClick;...
                        [fn*ones(numClicks,1),allClickOverallLickRate,...
                        period*ones(numClicks,1),MissOrHit]];
                    % store the click times
                    allClickTimes=[allClickTimes;[fn*ones(numClicks,1),...
                        clickArray,period*ones(numClicks,1),MissOrHit]];
                    %%% wheel speed
                    if szTiming(szCount,:).hasWheel
                        clickWheelData = ProcessWheelData_220630...
                            (clickArray,wheelData,params);
                        % smooth aburpt changes,
                        % 'sgolay':Savitzky-Golay filter
                        % smooths according to a quadratic polynomial that  
                        % is fitted over each window of data. 
                        % This method can be more effective than other 
                        % methods when the data varies rapidly.
                        clickWheelData_smoothed=cellfun(@(v)smoothdata...
                            (v,'sgolay'),clickWheelData,...
                            'UniformOutput',false);
                        % calculate the speed
                        % John D'Errico (2022). Movingslope 
                        % (https://www.mathworks.com/matlabcentral/fileexchange/16997-movingslope)
                        % MATLAB Central File Exchange. 
                        % Retrieved 06/30,2022.
                        clickWheelSpeed=...
                            cellfun(@(v)abs(movingslope(v,5,2)),...
                            clickWheelData_smoothed,'UniformOutput',false);
                        % Remove outliers
                        % (maybe not necessary if we smooth?)
                        % Fills with the nearest non-outlier value
                        clickWheelSpeed=cellfun(@(v)filloutliers...
                            (v,'nearest'),clickWheelSpeed,...
                            'UniformOutput',false);
                        WheelSpeed.AllClicks = ...
                            [WheelSpeed.AllClicks;clickWheelSpeed];
                        clickWheelSpeedAvg = ...
                            cellfun(@(v)mean(v,'omitnan'),clickWheelSpeed);
                        clickWheelSpeedMax = ...
                            cellfun(@(v)max(v,[],'omitnan'),...
                            clickWheelSpeed);             
                    else
                        WheelSpeed.AllClicks = [WheelSpeed.AllClicks;...
                            num2cell(nan(numClicks,1))];
                        [clickWheelSpeedAvg, clickWheelSpeedMax] = ...
                            deal(nan(numClicks,1));
                    end
                    WheelSpeed.Max.PerClick = [WheelSpeed.Max.PerClick;...
                        [fn*ones(numClicks,1),clickWheelSpeedMax,...
                        period*ones(numClicks,1),MissOrHit]];
                    WheelSpeed.Mean.PerClick = ...
                        [WheelSpeed.Mean.PerClick;...
                        [fn*ones(numClicks,1),clickWheelSpeedAvg,...
                        period*ones(numClicks,1),MissOrHit]];
                    % the time from sz onset for the seizure priod [array]
                    % the 1st column is time from seizure onset
                    % the 2nd is whether a click is a hit(1) or a miss(0)
                    if period == 2  % if this is a seizure period       
                        % calculate the time from sz onset for all clicks 
                        % and append to the All array
                        timeFromSzStart = [fn*ones(numClicks,1),...
                            (clickArray*params.SampleRate...
                            -(AllTimes(3)))./params.SampleRate,MissOrHit];
                        TimeFromSzStartAll = ...
                            [TimeFromSzStartAll;timeFromSzStart];
                    end
                    % number of hits and misses
                    hitCount = sum(hitInds);
                    missCount = sum(missInds);
                end
                % append to the global arrays that hold data for all files
                FirstLickAll.Miss{fn, period} = missFirstLick;
                FirstLickAll.Hit{fn, period} = hitFirstLick;
                FirstLickAll.AllClicks = [FirstLickAll.AllClicks;...
                    [fn*ones(numClicks,1),responseTimes,...
                    period*ones(numClicks,1),MissOrHit]];
                %
                NumLicksAll.Miss{fn, period} = missNumLicks;
                NumLicksAll.Hit{fn, period} = hitNumLicks;
                NumLicksAll.AllClicks = [NumLicksAll.AllClicks;...
                    [fn*ones(numClicks,1),allClickNumLicks,...
                    period*ones(numClicks,1),MissOrHit]];
                %
                ClickResponses(szCount,1) = fn;
                ClickResponses(szCount,(2:6)+(5*(period-1))) = ...
                    [mean(hitFirstLick,'omitnan'),hitCount,...
                    mean(missFirstLick,'omitnan'),missCount,...
                    (hitCount/(hitCount+missCount))];
                % 
                TimeCourse.Hit(:,szCount,period) = ...
                    mean(hitTimeCourses,1,'omitnan')';
                TimeCourse.Miss(:,szCount,period) = ...
                    mean(missTimeCourses,1,'omitnan')';
                TimeCourse.All(:,szCount,period) = ...
                    mean([hitTimeCourses;missTimeCourses],'omitnan');
                %
                LickRate.Max.PerSeizure(szCount,period) = ...
                    max(TimeCourse.All(binTimePts>=maxLickRateWindow(1)...
                    &binTimePts<=maxLickRateWindow(2),...
                    szCount,period))/binSize;
                %
                LickRate.Overall.PerSeizure(szCount,period) = ...
                    mean(LickRate.Overall.PerClick...
                    ((LickRate.Overall.PerClick(:,1)==szCount)&...
                    (LickRate.Overall.PerClick(:,3)==period),2));
                %
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
        Behavior.ByClick.AnimalsAll = Animal(allClickTimes(:,1));
    else
        Behavior.ByClick.AnimalsAll = [];
    end
    Behavior.ByClick.WheelSpeed.Max = WheelSpeed.Max.PerClick;
    Behavior.ByClick.WheelSpeed.Mean = WheelSpeed.Mean.PerClick;
    Behavior.ByClick.WheelSpeed.All = WheelSpeed.AllClicks;
    %%% Averages by seizures
    % Array with 4 colomns 
    % (1) sz file index, 
    % (2) number of clicks in the sz period of that file, 
    % (3) average time from sz onset across clicks for that sz,
    % (4) hit rate (hit clicks percentage) for the sz period in that file
    Behavior.BySeizure.TimeFromSzStart = GetAvgByClass_220704...
        (TimeFromSzStartAll(:,2:end),TimeFromSzStartAll(:,1));
    Behavior.BySeizure.IncludedSzDuration = IncludedSzDuration;
    Behavior.BySeizure.TimeCourse = TimeCourse;
    % 
    Behavior.BySeizure.MeanTimeCourse.Hit = ...
        [mean(TimeCourse.Hit(:,:,1),2,'omitnan'),...
        mean(TimeCourse.Hit(:,:,2),2,'omitnan'),...
        mean(TimeCourse.Hit(:,:,3),2,'omitnan'),...
        mean(TimeCourse.Hit(:,:,4),2,'omitnan')];
    Behavior.BySeizure.MeanTimeCourse.Miss = ...
        [mean(TimeCourse.Miss(:,:,1),2,'omitnan'),...
        mean(TimeCourse.Miss(:,:,2),2,'omitnan'), ...
        mean(TimeCourse.Miss(:,:,3),2,'omitnan'),...
        mean(TimeCourse.Miss(:,:,4),2,'omitnan')];
    %
    Behavior.BySeizure.HitRate = ClickResponses(:, 6:5:21);
    Behavior.BySeizure.HitLatency = ClickResponses(:,2:5:17);
    Behavior.BySeizure.MissLatency = ClickResponses(:,4:5:19);
    Behavior.BySeizure.DurationAndHitRate = ...
        [ClickResponses(:,11),IncludedSzDuration,ClickResponses(:,11)];
    %
    Behavior.BySeizure.ClickIntervals = avgClickIntervals;
    Behavior.BySeizure.NumClicks = NumClicksAll;
    Behavior.BySeizure.LickRateAll.Max = LickRate.Max.PerSeizure;
    Behavior.BySeizure.LickRateAll.Overall = LickRate.Overall.PerSeizure;
    %%% Averages by animal
    Behavior.ByAnimal.IncludedSzDuration = GetAvgByClass_220704...
        (IncludedSzDuration(:,2),Animal(IncludedSzDuration(:,1)));
    Behavior.ByAnimal.HitRate = GetAvgByClass_220704...
        (ClickResponses(:,6:5:21),Animal(seizures'));
    Behavior.ByAnimal.HitLatency = GetAvgByClass_220704...
        (ClickResponses(:,2:5:17),Animal(seizures'));
    Behavior.ByAnimal.DurationAndHitRate = GetAvgByClass_220704...
        ([IncludedSzDuration(:,2),...
        ClickResponses(:,11)],Animal(ClickResponses(:,1)));
    % Time course's dimension: (2+nBins)*nAnimals
    % The first two row: Animal and click numbers associated with
    %                    that animal in a certain period.
    % Row 2:(nBins+2) correspond to the actual lick rate time course. 
    % Each column corresponds to one animal.
    % Hit
    Behavior.ByAnimal.TimeCourse.Hit.Baseline = GetAvgByClass_220704...
        (TimeCourse.Hit(:,:,1)',Animal(seizures)')';
    Behavior.ByAnimal.TimeCourse.Hit.Ictal = GetAvgByClass_220704...
        (TimeCourse.Hit(:,:,2)',Animal(seizures)')';
    Behavior.ByAnimal.TimeCourse.Hit.PostIctal = GetAvgByClass_220704...
        (TimeCourse.Hit(:,:,3)',Animal(seizures)')';
    Behavior.ByAnimal.TimeCourse.Hit.Recovery = GetAvgByClass_220704...
        (TimeCourse.Hit(:,:,4)',Animal(seizures)')';
    % Miss
    Behavior.ByAnimal.TimeCourse.Miss.Baseline = GetAvgByClass_220704...
        (TimeCourse.Miss(:,:,1)',Animal(seizures)')';
    Behavior.ByAnimal.TimeCourse.Miss.Ictal = GetAvgByClass_220704...
        (TimeCourse.Miss(:,:,2)',Animal(seizures)')';
    Behavior.ByAnimal.TimeCourse.Miss.PostIctal = GetAvgByClass_220704...
        (TimeCourse.Miss(:,:,3)',Animal(seizures)')';
    Behavior.ByAnimal.TimeCourse.Miss.Recovery = GetAvgByClass_220704...
        (TimeCourse.Miss(:,:,4)',Animal(seizures)')';
    %
    Behavior.ByAnimal.NumClicks = GetAvgByClass_220704...
        (NumClicksAll,Animal(seizures));
    % save
    save(fullfile(SaveDir,'Behavior.mat'),'Behavior');
    fprintf('Done!\n'); 
end
