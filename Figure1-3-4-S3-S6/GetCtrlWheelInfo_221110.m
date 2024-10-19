% Get period based control wheel slope data
% Edit by Jerry, 11/10/2022
function [wheelSopeData,CtrlInclude,CtrlExclude,wheelSopeInclude] = ...
    GetCtrlWheelInfo_221110(params,szCell,szDir,szTiming,BLthd,SaveDir)
    sps = params.SampleRate; % 1000
    PreCtrlOffset = params.PreIctalOffset*sps; % baseline 60s*1000
    PostCtrlOffset = params.PostIctalOffset*sps; % postCtrl 30s*1000
    RecoveryPeriod = params.RecoveryPeriod*sps; % recovery 60s*1000
    % include/exclude array
    toInclude = string(szCell(2:end,3));
    szCount = 0; % counter for the number of included .mat data
    wheelCount = 0;
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
        if szTiming(szCount,:).hasWheel
            if szTiming(szCount,:).hasControl
                wheelCount = wheelCount+1;
                wAllTimes(wheelCount,:) = szTiming(szCount,:).AllTimes;
                wheelData(wheelCount,:) = szTiming(szCount,:).Wheel; 
            end
        end
    end
    %% wheel data should start before control baseline onset and ends after 
    for i = 1:wheelCount
        wheelTemp = wheelData(i,:);
        xx(i,:) = 0<(wAllTimes(i,9)-PreCtrlOffset-wheelTemp.start);
        yy(i,:) = wheelTemp.length>(wAllTimes(i,10)...
            +PostCtrlOffset+RecoveryPeriod-wheelTemp.start);
    end
    %% exclude xx=0 or yy = 0
    zz1 = find(xx==0);     zz2 = find(yy==0);
    zz = unique(cat(1,zz1,zz2));
    wAllTimes(zz,:)=[];     wheelData(zz,:)=[];
    %% Get period-based wheel slope data
    CtrlInclude = [];
    CtrlExclude = [];
    countInclude = 0;
    countExclude = 0;
    for i = 1:length(wheelData)
        disp(i);
        wheelTemp = wheelData(i,:);
        wheelvalue = wheelTemp.values;
        % find (0 <-> 5V) transition points 
        dwhv = diff(wheelvalue);
        dwhind = find(dwhv>1000);
        % make transitions continuous 
        for j = 1:length(dwhind)-1
            wheelvalue(dwhind(j):dwhind(j+1)) = ...
                wheelvalue(dwhind(j):dwhind(j+1))-5000*j;
        end
        % exclusion using max mice speed: 60m/min = 6000cm/60s = 100cm/s
        dwhv2 = diff(wheelvalue);
        dwhind2 = find(dwhv2>300 | dwhv2<-300); % choose value > 100
        wheelvalue(dwhind2)=nan;
        % get wheel slope, abs make it no negative slope
        wheelSlopeTemp = abs(movingslope(wheelvalue,5,2));
        % get rid of peaks
        wheelSlope = filloutliers(wheelSlopeTemp,'nearest');
        % control baseline wheelSlope
        wBaselineTime = round(((wAllTimes(i,9)-PreCtrlOffset+1):...
            wAllTimes(i,9))-wheelTemp.start);
        wBaselineData(i,:) = wheelSlope(wBaselineTime(1:PreCtrlOffset));
        % run/still based on threshold for baseline
        wBLD_Temp = mean(wBaselineData(i,:),'omitnan');
        if wBLD_Temp>BLthd
            countInclude = countInclude+1;
            CtrlInclude(countInclude) = i;
        else
            countExclude = countExclude+1;
            CtrlExclude(countExclude) = i;
        end
        % Ctrl wheelSlope whose length equals to sz duration
        wCtrlTime = round(((wAllTimes(i,9)+1):wAllTimes(i,10))-...
            wheelTemp.start);
        wCtrlData{i} = wheelSlope(wCtrlTime);
        % postCtrl wheelSlope
        wPostCtrlTime = round(((wAllTimes(i,10)+1):...
            (wAllTimes(i,10)+PostCtrlOffset))-wheelTemp.start);
        wPostCtrlData(i,:)=wheelSlope(wPostCtrlTime(1:PostCtrlOffset));
        % recovery wheelSlope
        wRecoveryTime = round(((wAllTimes(i,10)+PostCtrlOffset+1)...
            :(wAllTimes(i,10)+PostCtrlOffset+RecoveryPeriod))-...
            wheelTemp.start);
        wRecoveryData(i,:) = wheelSlope(wRecoveryTime(1:RecoveryPeriod));
    end
    %% mean across time
    wBaseline_mean = mean(wBaselineData,2,'omitnan');
    for i = 1:length(wCtrlData)
        wCtrl_mean(i,:) = mean(wCtrlData{i},'omitnan');
    end
    wPostCtrl_mean = mean(wPostCtrlData,2,'omitnan');
    wRecovery_mean = mean(wRecoveryData,2,'omitnan');
    wheelSopeData = ...
        cat(2,wBaseline_mean,wCtrl_mean,wPostCtrl_mean,wRecovery_mean);
    wheelSopeInclude = wheelSopeData(CtrlInclude,:);
    % save .mat file
    save(fullfile(SaveDir,'wheelCtrlSopeData.mat'),'wheelSopeData'); 
    save(fullfile(SaveDir,'WheelCtrlInclude.mat'),'CtrlInclude'); 
    save(fullfile(SaveDir,'WheelCtrlExclude.mat'),'CtrlExclude'); 
    save(fullfile(SaveDir,'wheelCtrlSopeInclude.mat'),'wheelSopeInclude'); 
    disp('done!');
end