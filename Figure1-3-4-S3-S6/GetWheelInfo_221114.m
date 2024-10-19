% Get period based wheel slope data
% Edited by Jerry. 11/14/2022
function [wheelSopeData,szInclude,szExclude,wheelSopeInclude] = ...
    GetWheelInfo_221114(params,szCell,szDir,szTiming,BLthd,SaveDir)
    sps = params.SampleRate; % 1000
    PreIctalOffset = params.PreIctalOffset*sps; % baseline 60s*1000
    IctalCutoff = params.IctalCutoff*sps; % 15s*1000 ictal cutoff
    szTimeSkip = params.SeizureOnsetTimeSkip*sps; % 2.5s*1000
    PostIctalOffset = params.PostIctalOffset*sps; % postictal 30s*1000
    RecoveryPeriod = params.RecoveryPeriod*sps; % recovery 60s*1000
    % include/exclude array
    toInclude = string(szCell(2:end,3));
    szCount = 0; % counter for the number of included seizures
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
            wheelCount = wheelCount+1;
            wAllTimes(wheelCount,:) = szTiming(szCount,:).AllTimes;
            wheelData(wheelCount,:) = szTiming(szCount,:).Wheel; %[struct]
        end
    end
    %% check that wheel data start before baseline onset and ends after 
    for i = 1:wheelCount
        wheelTemp = wheelData(i,:);
        % 3s comes from the sz current train
        xx(i,:) = 0<(wAllTimes(i,1)-3000-wheelTemp.start);
        yy(i,:) = wheelTemp.length>(max(wAllTimes(i,:))-wheelTemp.start);
    end
    %% exclude xx=0 or yy = 0
    zz1 = find(xx==0);     zz2 = find(yy==0);
    zz = unique(cat(1,zz1,zz2));
    wAllTimes(zz,:)=[];     wheelData(zz,:)=[];
    %% Get period-based wheel slope data
    szInclude = [];
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
        % baseline wheelSlope, 3000 is used to get rid of 2s current train
        wBaselineTime = round...
            (((wAllTimes(i,1)-3000):(wAllTimes(i,2)-2999))-wheelTemp.start);
        wBaselineData(i,:) = wheelSlope(wBaselineTime(1:PreIctalOffset));
        %% remove 
        wBLD_Temp = mean(wBaselineData(i,:),'omitnan');
        if wBLD_Temp>BLthd
            countInclude = countInclude+1;
            szInclude(countInclude) = i;
        else
            countExclude = countExclude+1;
            szExclude(countExclude) = i;
        end
        % ictal wheelSlope
        wIctalTime = round(((wAllTimes(i,3)+szTimeSkip):...
            (wAllTimes(i,4)+1))-wheelTemp.start);
        wIctalData{i} = wheelSlope(wIctalTime);
        % postictal wheelSlope
        wPostictalTime = ...
            round((wAllTimes(i,5):(wAllTimes(i,6)+1))-wheelTemp.start);
        wPostictalData(i,:)=wheelSlope(wPostictalTime(1:PostIctalOffset));
        % recovery wheelSlope
        wRecoveryTime = ...
            round((wAllTimes(i,7):(wAllTimes(i,8)+1))-wheelTemp.start);
        wRecoveryData(i,:) = wheelSlope(wRecoveryTime(1:RecoveryPeriod));
    end
    %% mean across time
    wBaseline_mean = mean(wBaselineData,2,'omitnan');
    for i = 1:length(wIctalData)
        wIctal_mean(i,:) = mean(wIctalData{i},'omitnan');
    end 
    wPostictal_mean = mean(wPostictalData,2,'omitnan');
    wRecovery_mean = mean(wRecoveryData,2,'omitnan');
    wheelSopeData = ...
        cat(2,wBaseline_mean,wIctal_mean,wPostictal_mean,wRecovery_mean);
    wheelSopeInclude = wheelSopeData(szInclude,:);
    % save .mat file
    save(fullfile(SaveDir,'wheelSopeData.mat'),'wheelSopeData'); 
    save(fullfile(SaveDir,'WheelInclude.mat'),'szInclude'); 
    save(fullfile(SaveDir,'WheelExclude.mat'),'szExclude'); 
    save(fullfile(SaveDir,'wheelSopeInclude.mat'),'wheelSopeInclude'); 
    disp('done!');
end