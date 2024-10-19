% Click based spectrogram/psd/bandpower analysis
% Edited by Jerry. 10/12/2022
function [BClk_SPGM,BClk_PSD,BClk_BP,IClk_SPGM,IClk_PSD,IClk_BP,...
    PClk_SPGM,PClk_PSD,PClk_BP,RClk_SPGM,RClk_PSD,RClk_BP]=...
    AnalyzeClickFrequencyV1_221012(params,Behavior,szCell,szDir,szTiming,SaveDir)
    %% Get relevant parameters 
    % electrophysiology analysis pre and post stim times 
    PreStimOffset = params.EPPreStimOffset; 
    PostStimOffset = params.EPPostStimOffset; 
    sps = params.SampleRate; % sample rate 
    szTimeSkip = params.SeizureOnsetTimeSkip; % sz onset time skip [s]
    artifactPeriod = params.StimulationArtifactPeriod*sps; % artifact [ms]
    % file numbers 
    fileNums = sort(unique(Behavior.ByClick.ClickTimes(:, 1))); 
    % include/exclude array for all .mat files
    toInclude = string(szCell(2:end,3));
    % number of seizures included
    szCount = 0; 
    for fn = 1:length(szDir)
        filepath = szDir{fn}; 
        disp(filepath);
        if strcmp(toInclude(fn), 'Include')
            szCount = szCount+1;  % increment the number of seizures 
        else
            continue;   % move on to next iteration 
        end
        % indeces of baseline,ictal,postical,recovery,and control
        AllTimesInds = szTiming(szCount).AllTimesInds;
        % LFP length is the sz .mat file length
        LO_LFP = szTiming(szCount).LO_LFP; 
        HC_LFP = szTiming(szCount).HC_LFP; 
        hasControl = szTiming(szCount).hasControl;
        % remove artifact following sz/control stimulation using polyfit
        % LO LFP artifact remove, degree 5
        LO_LFP_artifact = ...
            LO_LFP.values(AllTimesInds(3):AllTimesInds(3)+artifactPeriod);
        LO_LFP.values(AllTimesInds(3):AllTimesInds(3)+artifactPeriod) = ...
            RemoveArtifacts_220802(LO_LFP_artifact,5);
        % HC LFP artifact remove, degree 7
        HC_LFP_artifact = ...
            HC_LFP.values(AllTimesInds(3):AllTimesInds(3)+artifactPeriod);
        HC_LFP.values(AllTimesInds(3):AllTimesInds(3)+artifactPeriod) = ...
            RemoveArtifacts_220802(HC_LFP_artifact,7); 
        % 1uA control stimulation LO/HC LFP artifact 
        if hasControl
            % control LO LFP artifact remove, degree 5
            ControlLO_LFP_artifact = LO_LFP.values...
                (AllTimesInds(9):AllTimesInds(9)+artifactPeriod); 
            LO_LFP.values...
                (AllTimesInds(9):AllTimesInds(9)+artifactPeriod) =...
                RemoveArtifacts_220802(ControlLO_LFP_artifact,5); 
            % control HC LFP artifact remove, degree 7
            ControlHC_LFP_artifact = HC_LFP.values...
                (AllTimesInds(9):AllTimesInds(9)+artifactPeriod); 
            HC_LFP.values...
                (AllTimesInds(9):AllTimesInds(9)+artifactPeriod) =...
                RemoveArtifacts_220802(ControlHC_LFP_artifact,7); 
        end
        % skip clicks that are within the seizure onset skip period
        % click times for the clicks in this file 
        fileClickInfo = Behavior.ByClick.ClickTimes...
            (Behavior.ByClick.ClickTimes(:,1)==fileNums(fn),:); 
        fileClickTimes = fileClickInfo(:,2);
        % seizure onset and end [s]
        szStart = szTiming(szCount).AllTimes(3)/sps; 
        % skip clicks that are within the seizure onset skip period
        skipCondition1 = (fileClickTimes>szStart) & ...
            ((fileClickTimes-PreStimOffset)<=(szStart+szTimeSkip)); 
        % skip clicks that don't have enough pre stimulus data  &&
        skipCondition2 = round((fileClickTimes-PreStimOffset).*sps)<1;
        skipCondition3 = round((fileClickTimes+PostStimOffset).*sps)>...
            length(szTiming(szCount).LO_LFP.values); 
        clicksToSkip = skipCondition1|skipCondition2|skipCondition3;
        % LO,HC LFP
        LO_values = LO_LFP.values;
        HC_values = HC_LFP.values;
        % get reference spectrogram,psd,and bandpower from baseline
        firstClickIctal = 1; 
        Refcn = 0;
        for cn = 1:length(fileClickTimes)
            % deal with first click in ictal period
            if ~clicksToSkip(cn) && sum(skipCondition1)==0
                if fileClickInfo(cn, 3)==2 && fileClickTimes(cn)>szStart
                    if firstClickIctal
                        clicksToSkip(cn) = 1;
                        firstClickIctal = 0; 
                    end
                end
            end                
            if ~clicksToSkip(cn)            
                % baseline click
                if fileClickInfo(cn,3)==1
                    Refcn = Refcn + 1;
                    BStart = round((fileClickTimes(cn)-PreStimOffset)*sps);
                    BEnd = round((fileClickTimes(cn)+PostStimOffset)*sps);
                    BClkLO = LO_values(BStart+1:BEnd);
                    BClkHC = HC_values(BStart+1:BEnd);
                    % Reference click LO
                    [RefClkfy,RefClktx,RefClkLO_SPGM(:,:,Refcn),...
                        RefClkfw,RefClkLO_PSD(:,Refcn),...
                        RefClkLO_deltaBP(:,Refcn),RefClkLO_betaBP(:,Refcn)] = ...
                        GetFrequencyAnalysis_221011(params,BClkLO);
                    % Reference click HC
                    [~,~,RefClkHC_SPGM(:,:,Refcn),~,RefClkHC_PSD(:,Refcn),...
                        RefClkHC_deltaBP(:,Refcn),RefClkHC_betaBP(:,Refcn)] = ...
                        GetFrequencyAnalysis_221011(params,BClkHC);
                end
            end
        end
        % get reference spectrogram for this sz.mat, linear units
        RefClkLO_SPGM_mean = mean(RefClkLO_SPGM,3,'omitnan');
        RefClkLO_SPGM_mean = mean(RefClkLO_SPGM_mean,2,'omitnan');
        RefClkHC_SPGM_mean = mean(RefClkHC_SPGM,3,'omitnan');
        RefClkHC_SPGM_mean = mean(RefClkHC_SPGM_mean,2,'omitnan');
        % get reference psd for this sz.mat, linear units
        RefClkLO_PSD_mean = mean(RefClkLO_PSD,2,'omitnan');
        RefClkHC_PSD_mean = mean(RefClkHC_PSD,2,'omitnan');
        % reference LO bandpower for this sz.mat, linear units
        RefLO_deltaBP_mean = mean(RefClkLO_deltaBP,2,'omitnan'); % delta
        RefLO_betaBP_mean = mean(RefClkLO_betaBP,2,'omitnan'); % beta
        % reference HC bandpower for this sz.mat, linear units
        RefHC_deltaBP_mean = mean(RefClkHC_deltaBP,2,'omitnan'); % delta
        RefHC_betaBP_mean = mean(RefClkHC_betaBP,2,'omitnan'); % beta
        %
        % get periods and normed periods spectrogram,psd,and bandpower
        firstClickIctal = 1; 
        % Baseline
        HBcn = 0; % Baseline counter, all hits
        %
        xl = length(RefClkfy);
        yl = PreStimOffset+PostStimOffset;
        % hit ictal LO
        HIcn = 0; % Hit Ictal counter
        NormHIClkLO_SPGM = nan(xl,yl); % Normed Hit ictal LO spectrogram
        NormHIClkLO_PSD = nan(xl,1); % Normed Hit ictal LO psd
        NormHIClkLO_deltaBP = nan; % Normed Hit ictal LO deltaBP
        NormHIClkLO_betaBP = nan; % Normed Hit ictal LO betaBP
        % hit ictal HC
        NormHIClkHC_SPGM = nan(xl,yl); % Normed Hit ictal HC spectrogram
        NormHIClkHC_PSD = nan(xl,1); % Normed Hit ictal HC psd
        NormHIClkHC_deltaBP = nan; % Normed Hit ictal HC deltaBP
        NormHIClkHC_betaBP = nan; % Normed Hit ictal HC betaBP
        % miss ictal LO
        MIcn = 0; % Miss Ictal counter
        NormMIClkLO_SPGM = nan(xl,yl); % Normed Miss ictal LO spectrogram
        NormMIClkLO_PSD = nan(xl,1); % Normed Hit ictal LO psd
        NormMIClkLO_deltaBP = nan; % Normed Hit ictal LO deltaBP
        NormMIClkLO_betaBP = nan; % Normed Hit ictal LO betaBP        
        % miss ictal HC
        NormMIClkHC_SPGM = nan(xl,yl); % Normed Miss ictal HC spectrogram
        NormMIClkHC_PSD = nan(xl,1); % Normed Hit ictal HC psd
        NormMIClkHC_deltaBP = nan; % Normed Hit ictal HC deltaBP
        NormMIClkHC_betaBP = nan; % Normed Hit ictal HC betaBP
        %
        % hit Postictal LO
        HPcn = 0; % hit Postictal counter
        NormHPClkLO_SPGM = nan(xl,yl); % Normed Hit Postictal LO spectrogram
        NormHPClkLO_PSD = nan(xl,1); % Normed Hit Postictal LO psd
        NormHPClkLO_deltaBP = nan; % Normed Hit Postictal LO deltaBP
        NormHPClkLO_betaBP = nan; % Normed Hit Postictal LO betaBP
        % hit Postictal HC
        NormHPClkHC_SPGM = nan(xl,yl); % Normed Hit Postictal HC spectrogram
        NormHPClkHC_PSD = nan(xl,1); % Normed Hit Postictal HC psd
        NormHPClkHC_deltaBP = nan; % Normed Hit Postictal HC deltaBP
        NormHPClkHC_betaBP = nan; % Normed Hit Postictal HC betaBP
        % miss Postictal LO
        MPcn = 0; % miss Postictal counter
        NormMPClkLO_SPGM = nan(xl,yl); % Normed miss Postictal LO spectrogram
        NormMPClkLO_PSD = nan(xl,1); % Normed miss Postictal LO psd
        NormMPClkLO_deltaBP = nan; % Normed miss Postictal LO deltaBP
        NormMPClkLO_betaBP = nan; % Normed miss Postictal LO betaBP
        % miss Postictal HC
        NormMPClkHC_SPGM = nan(xl,yl); % Normed miss Postictal HC spectrogram
        NormMPClkHC_PSD = nan(xl,1); % Normed miss Postictal HC psd
        NormMPClkHC_deltaBP = nan; % Normed miss Postictal HC deltaBP
        NormMPClkHC_betaBP = nan; % Normed miss Postictal HC betaBP
        %
        % hit recovery LO
        HRcn = 0; % hit recovery counter
        NormHRClkLO_SPGM = nan(xl,yl); % Normed Hit recovery LO spectrogram
        NormHRClkLO_PSD = nan(xl,1); % Normed Hit recovery LO psd
        NormHRClkLO_deltaBP = nan; % Normed Hit recovery LO deltaBP
        NormHRClkLO_betaBP = nan; % Normed Hit recovery LO betaBP
        % hit recovery HC
        NormHRClkHC_SPGM = nan(xl,yl); % Normed Hit recovery HC spectrogram
        NormHRClkHC_PSD = nan(xl,1); % Normed Hit recovery HC psd
        NormHRClkHC_deltaBP = nan; % Normed Hit recovery HC deltaBP
        NormHRClkHC_betaBP = nan; % Normed Hit recovery HC betaBP
        % miss recovery LO
        MRcn = 0; % miss recovery counter
        NormMRClkLO_SPGM = nan(xl,yl); % Normed miss recovery LO spectrogram
        NormMRClkLO_PSD = nan(xl,1); % Normed miss recovery LO psd
        NormMRClkLO_deltaBP = nan; % Normed miss recovery LO deltaBP
        NormMRClkLO_betaBP = nan; % Normed miss recovery LO betaBP
        % miss recovery HC
        NormMRClkHC_SPGM = nan(xl,yl); % Normed miss recovery HC spectrogram
        NormMRClkHC_PSD = nan(xl,1); % Normed miss recovery HC psd
        NormMRClkHC_deltaBP = nan; % Normed miss recovery HC deltaBP
        NormMRClkHC_betaBP = nan; % Normed miss recovery HC betaBP
        % look clicks
        for cn = 1:length(fileClickTimes)
            % deal with first click in ictal period
            if ~clicksToSkip(cn) && sum(skipCondition1)==0
                if fileClickInfo(cn, 3)==2 && fileClickTimes(cn)>szStart
                    if firstClickIctal
                        clicksToSkip(cn) = 1;
                        firstClickIctal = 0; 
                    end
                end
            end                
            if ~clicksToSkip(cn)
                % baseline click, all hits
                if fileClickInfo(cn,3)==1
                    HBcn = HBcn + 1;
                    BStart = round((fileClickTimes(cn)-PreStimOffset)*sps);
                    BEnd = round((fileClickTimes(cn)+PostStimOffset)*sps);
                    BClkLO = LO_values(BStart+1:BEnd);
                    BClkHC = HC_values(BStart+1:BEnd);
                    % baseline click LO
                    [BClkfy,BClktx,BClkLO_SPGM,BClkfw,BClkLO_PSD,...
                        BClkLO_deltaBP,BClkLO_betaBP] = ...
                        GetFrequencyAnalysis_221011(params,BClkLO);
                    % Normalization
                    NormBClkLO_SPGM(:,:,HBcn) = BClkLO_SPGM./RefClkLO_SPGM_mean;
                    NormBClkLO_PSD(:,HBcn) = BClkLO_PSD./RefClkLO_PSD_mean;
                    NormBClkLO_deltaBP(:,HBcn) = BClkLO_deltaBP/RefLO_deltaBP_mean;
                    NormBClkLO_betaBP(:,HBcn) = BClkLO_betaBP/RefLO_betaBP_mean;
                    %
                    % baseline click HC
                    [~,~,BClkHC_SPGM,~,BClkHC_PSD,BClkHC_deltaBP,BClkHC_betaBP] = ...
                        GetFrequencyAnalysis_221011(params,BClkHC);
                    % Normalization
                    NormBClkHC_SPGM(:,:,HBcn) = BClkHC_SPGM./RefClkHC_SPGM_mean;
                    NormBClkHC_PSD(:,HBcn) = BClkHC_PSD./RefClkHC_PSD_mean;
                    NormBClkHC_deltaBP(:,HBcn) = BClkHC_deltaBP/RefHC_deltaBP_mean;
                    NormBClkHC_betaBP(:,HBcn) = BClkHC_betaBP/RefHC_betaBP_mean;
                end
                % ictal period
                if fileClickInfo(cn,3)==2
                    % hit clicks
                    if fileClickInfo(cn,4) == 1
                        HIcn = HIcn + 1;
                        HIStart = round((fileClickTimes(cn)-PreStimOffset)*sps);
                        HIEnd = round((fileClickTimes(cn)+PostStimOffset)*sps);
                        HIClkLO = LO_values(HIStart+1:HIEnd);
                        HIClkHC = HC_values(HIStart+1:HIEnd);
                        % hit ictal click LO
                        [IClkfy,IClktx,HIClkLO_SPGM,IClkfw,HIClkLO_PSD,...
                            HIClkLO_deltaBP,HIClkLO_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,HIClkLO);
                        % Normalization
                        NormHIClkLO_SPGM(:,:,HIcn) = HIClkLO_SPGM./RefClkLO_SPGM_mean;
                        NormHIClkLO_PSD(:,HIcn) = HIClkLO_PSD./RefClkLO_PSD_mean;
                        NormHIClkLO_deltaBP(:,HIcn) = HIClkLO_deltaBP/RefLO_deltaBP_mean;
                        NormHIClkLO_betaBP(:,HIcn) = HIClkLO_betaBP/RefLO_betaBP_mean;
                        % hit ictal click HC
                        [~,~,HIClkHC_SPGM,~,HIClkHC_PSD,HIClkHC_deltaBP,HIClkHC_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,HIClkHC);
                        % Normalization
                        NormHIClkHC_SPGM(:,:,HIcn) = HIClkHC_SPGM./RefClkHC_SPGM_mean;
                        NormHIClkHC_PSD(:,HIcn) = HIClkHC_PSD./RefClkHC_PSD_mean;
                        NormHIClkHC_deltaBP(:,HIcn) = HIClkHC_deltaBP/RefHC_deltaBP_mean;
                        NormHIClkHC_betaBP(:,HIcn) = HIClkHC_betaBP/RefHC_betaBP_mean;
                    % miss clicks
                    else 
                        MIcn = MIcn + 1;
                        MIStart = round((fileClickTimes(cn)-PreStimOffset)*sps);
                        MIEnd = round((fileClickTimes(cn)+PostStimOffset)*sps);
                        MIClkLO = LO_values(MIStart+1:MIEnd);
                        MIClkHC = HC_values(MIStart+1:MIEnd);
                        % miss ictal click LO
                        [IClkfy,IClktx,MIClkLO_SPGM,IClkfw,MIClkLO_PSD,...
                            MIClkLO_deltaBP,MIClkLO_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,MIClkLO);
                        % Normalization
                        NormMIClkLO_SPGM(:,:,MIcn) = MIClkLO_SPGM./RefClkLO_SPGM_mean;
                        NormMIClkLO_PSD(:,MIcn) = MIClkLO_PSD./RefClkLO_PSD_mean;
                        NormMIClkLO_deltaBP(:,MIcn) = MIClkLO_deltaBP/RefLO_deltaBP_mean;
                        NormMIClkLO_betaBP(:,MIcn) = MIClkLO_betaBP/RefLO_betaBP_mean;
                        % miss ictal click HC
                        [~,~,MIClkHC_SPGM,~,MIClkHC_PSD,MIClkHC_deltaBP,MIClkHC_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,MIClkHC);
                        % Normalization
                        NormMIClkHC_SPGM(:,:,MIcn) = MIClkHC_SPGM./RefClkHC_SPGM_mean;
                        NormMIClkHC_PSD(:,MIcn) = MIClkHC_PSD./RefClkHC_PSD_mean;
                        NormMIClkHC_deltaBP(:,MIcn) = MIClkHC_deltaBP/RefHC_deltaBP_mean;
                        NormMIClkHC_betaBP(:,MIcn) = MIClkHC_betaBP/RefHC_betaBP_mean;
                    end
                end
                %
                % postictal period
                if fileClickInfo(cn,3)==3
                    % hit clicks
                    if fileClickInfo(cn,4) == 1
                        HPcn = HPcn + 1;
                        HPStart = round((fileClickTimes(cn)-PreStimOffset)*sps);
                        HPEnd = round((fileClickTimes(cn)+PostStimOffset)*sps);
                        HPClkLO = LO_values(HPStart+1:HPEnd);
                        HPClkHC = HC_values(HPStart+1:HPEnd);
                        % hit postictal click LO
                        [PClkfy,PClktx,HPClkLO_SPGM,PClkfw,HPClkLO_PSD,...
                            HPClkLO_deltaBP,HPClkLO_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,HPClkLO);
                        % Normalization
                        NormHPClkLO_SPGM(:,:,HPcn) = HPClkLO_SPGM./RefClkLO_SPGM_mean;
                        NormHPClkLO_PSD(:,HPcn) = HPClkLO_PSD./RefClkLO_PSD_mean;
                        NormHPClkLO_deltaBP(:,HPcn) = HPClkLO_deltaBP/RefLO_deltaBP_mean;
                        NormHPClkLO_betaBP(:,HPcn) = HPClkLO_betaBP/RefLO_betaBP_mean;
                        % hit postictal click HC
                        [~,~,HPClkHC_SPGM,~,HPClkHC_PSD,HPClkHC_deltaBP,HPClkHC_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,HPClkHC);
                        % Normalization
                        NormHPClkHC_SPGM(:,:,HPcn) = HPClkHC_SPGM./RefClkHC_SPGM_mean;
                        NormHPClkHC_PSD(:,HPcn) = HPClkHC_PSD./RefClkHC_PSD_mean;
                        NormHPClkHC_deltaBP(:,HPcn) = HPClkHC_deltaBP/RefHC_deltaBP_mean;
                        NormHPClkHC_betaBP(:,HPcn) = HPClkHC_betaBP/RefHC_betaBP_mean;
                    % miss clicks
                    else 
                        MPcn = MPcn + 1;
                        MPStart = round((fileClickTimes(cn)-PreStimOffset)*sps);
                        MPEnd = round((fileClickTimes(cn)+PostStimOffset)*sps);
                        MPClkLO = LO_values(MPStart+1:MPEnd);
                        MPClkHC = HC_values(MPStart+1:MPEnd);
                        % miss postictal click LO
                        [PClkfy,PClktx,MPClkLO_SPGM,PClkfw,MPClkLO_PSD,...
                            MPClkLO_deltaBP,MPClkLO_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,MPClkLO);
                        % Normalization
                        NormMPClkLO_SPGM(:,:,MPcn) = MPClkLO_SPGM./RefClkLO_SPGM_mean;
                        NormMPClkLO_PSD(:,MPcn) = MPClkLO_PSD./RefClkLO_PSD_mean;
                        NormMPClkLO_deltaBP(:,MPcn) = MPClkLO_deltaBP/RefLO_deltaBP_mean;
                        NormMPClkLO_betaBP(:,MPcn) = MPClkLO_betaBP/RefLO_betaBP_mean;
                        % miss postictal click HC
                        [~,~,MPClkHC_SPGM,~,MPClkHC_PSD,MPClkHC_deltaBP,MPClkHC_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,MPClkHC);
                        % Normalization
                        NormMPClkHC_SPGM(:,:,MPcn) = MPClkHC_SPGM./RefClkHC_SPGM_mean;
                        NormMPClkHC_PSD(:,MPcn) = MPClkHC_PSD./RefClkHC_PSD_mean;
                        NormMPClkHC_deltaBP(:,MPcn) = MPClkHC_deltaBP/RefHC_deltaBP_mean;
                        NormMPClkHC_betaBP(:,MPcn) = MPClkHC_betaBP/RefHC_betaBP_mean;
                    end
                end
                %
                % recovery period
                if fileClickInfo(cn,3)==4
                    % hit clicks
                    if fileClickInfo(cn,4) == 1
                        HRcn = HRcn + 1;
                        HRStart = round((fileClickTimes(cn)-PreStimOffset)*sps);
                        HREnd = round((fileClickTimes(cn)+PostStimOffset)*sps);
                        HRClkLO = LO_values(HRStart+1:HREnd);
                        HRClkHC = HC_values(HRStart+1:HREnd);
                        % hit recovery click LO
                        [RClkfy,RClktx,HRClkLO_SPGM,RClkfw,HRClkLO_PSD,...
                            HRClkLO_deltaBP,HRClkLO_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,HRClkLO);
                        % Normalization
                        NormHRClkLO_SPGM(:,:,HRcn) = HRClkLO_SPGM./RefClkLO_SPGM_mean;
                        NormHRClkLO_PSD(:,HRcn) = HRClkLO_PSD./RefClkLO_PSD_mean;
                        NormHRClkLO_deltaBP(:,HRcn) = HRClkLO_deltaBP/RefLO_deltaBP_mean;
                        NormHRClkLO_betaBP(:,HRcn) = HRClkLO_betaBP/RefLO_betaBP_mean;
                        % hit recovery click HC
                        [~,~,HRClkHC_SPGM,~,HRClkHC_PSD,HRClkHC_deltaBP,HRClkHC_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,HRClkHC);
                        % Normalization
                        NormHRClkHC_SPGM(:,:,HRcn) = HRClkHC_SPGM./RefClkHC_SPGM_mean;
                        NormHRClkHC_PSD(:,HRcn) = HRClkHC_PSD./RefClkHC_PSD_mean;
                        NormHRClkHC_deltaBP(:,HRcn) = HRClkHC_deltaBP/RefHC_deltaBP_mean;
                        NormHRClkHC_betaBP(:,HRcn) = HRClkHC_betaBP/RefHC_betaBP_mean;
                    % miss clicks
                    else 
                        MRcn = MRcn + 1;
                        MRStart = round((fileClickTimes(cn)-PreStimOffset)*sps);
                        MREnd = round((fileClickTimes(cn)+PostStimOffset)*sps);
                        MRClkLO = LO_values(MRStart+1:MREnd);
                        MRClkHC = HC_values(MRStart+1:MREnd);
                        % miss recovery click LO
                        [RClkfy,RClktx,MRClkLO_SPGM,RClkfw,MRClkLO_PSD,...
                            MRClkLO_deltaBP,MRClkLO_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,MRClkLO);
                        % Normalization
                        NormMRClkLO_SPGM(:,:,MRcn) = MRClkLO_SPGM./RefClkLO_SPGM_mean;
                        NormMRClkLO_PSD(:,MRcn) = MRClkLO_PSD./RefClkLO_PSD_mean;
                        NormMRClkLO_deltaBP(:,MRcn) = MRClkLO_deltaBP/RefLO_deltaBP_mean;
                        NormMRClkLO_betaBP(:,MRcn) = MRClkLO_betaBP/RefLO_betaBP_mean;
                        % miss recovery click HC
                        [~,~,MRClkHC_SPGM,~,MRClkHC_PSD,MRClkHC_deltaBP,MRClkHC_betaBP] = ...
                            GetFrequencyAnalysis_221011(params,MRClkHC);
                        % Normalization
                        NormMRClkHC_SPGM(:,:,MRcn) = MRClkHC_SPGM./RefClkHC_SPGM_mean;
                        NormMRClkHC_PSD(:,MRcn) = MRClkHC_PSD./RefClkHC_PSD_mean;
                        NormMRClkHC_deltaBP(:,MRcn) = MRClkHC_deltaBP/RefHC_deltaBP_mean;
                        NormMRClkHC_betaBP(:,MRcn) = MRClkHC_betaBP/RefHC_betaBP_mean;
                    end
                end
            end
        end
        % get normed baseline spectrogram for this sz.mat, linear units
        NBClkLO_SPGM_mat(:,:,szCount) = mean(NormBClkLO_SPGM,3,'omitnan');
        NBClkHC_SPGM_mat(:,:,szCount) = mean(NormBClkHC_SPGM,3,'omitnan');
        % get normed baseline psd for this sz.mat, linear units
        NBClkLO_PSD_mat(:,szCount) = mean(NormBClkLO_PSD,2,'omitnan');
        NBClkHC_PSD_mat(:,szCount) = mean(NormBClkHC_PSD,2,'omitnan');
        % normed baseline LO bandpower for this sz.mat, linear units
        NBClkLO_deltaBP_mat(:,szCount) = mean(NormBClkLO_deltaBP,2,'omitnan'); % delta
        NBClkLO_betaBP_mat(:,szCount) = mean(NormBClkLO_betaBP,2,'omitnan'); % beta
        % normed baseline HC bandpower for this sz.mat, linear units
        NBClkHC_deltaBP_mat(:,szCount) = mean(NormBClkHC_deltaBP,2,'omitnan'); % delta
        NBClkHC_betaBP_mat(:,szCount) = mean(NormBClkHC_betaBP,2,'omitnan'); % beta
        %
        % get normed Hit ictal spectrogram for this sz.mat, linear units
        NHIClkLO_SPGM_mat(:,:,szCount) = mean(NormHIClkLO_SPGM,3,'omitnan');
        NHIClkHC_SPGM_mat(:,:,szCount) = mean(NormHIClkHC_SPGM,3,'omitnan');
        % get normed Hit ictal psd for this sz.mat, linear units
        NHIClkLO_PSD_mat(:,szCount) = mean(NormHIClkLO_PSD,2,'omitnan');
        NHIClkHC_PSD_mat(:,szCount) = mean(NormHIClkHC_PSD,2,'omitnan');
        % normed Hit ictal LO bandpower for this sz.mat, linear units
        NHIClkLO_deltaBP_mat(:,szCount) = mean(NormHIClkLO_deltaBP,2,'omitnan'); % delta
        NHIClkLO_betaBP_mat(:,szCount) = mean(NormHIClkLO_betaBP,2,'omitnan'); % beta
        % normed Hit ictal HC bandpower for this sz.mat, linear units
        NHIClkHC_deltaBP_mat(:,szCount) = mean(NormHIClkHC_deltaBP,2,'omitnan'); % delta
        NHIClkHC_betaBP_mat(:,szCount) = mean(NormHIClkHC_betaBP,2,'omitnan'); % beta
        %
        % get normed Miss ictal spectrogram for this sz.mat, linear units
        NMIClkLO_SPGM_mat(:,:,szCount) = mean(NormMIClkLO_SPGM,3,'omitnan');
        NMIClkHC_SPGM_mat(:,:,szCount) = mean(NormMIClkHC_SPGM,3,'omitnan');
        % get normed Miss ictal psd for this sz.mat, linear units
        NMIClkLO_PSD_mat(:,szCount) = mean(NormMIClkLO_PSD,2,'omitnan');
        NMIClkHC_PSD_mat(:,szCount) = mean(NormMIClkHC_PSD,2,'omitnan');
        % normed Miss ictal LO bandpower for this sz.mat, linear units
        NMIClkLO_deltaBP_mat(:,szCount) = mean(NormMIClkLO_deltaBP,2,'omitnan'); % delta
        NMIClkLO_betaBP_mat(:,szCount) = mean(NormMIClkLO_betaBP,2,'omitnan'); % beta
        % normed Miss ictal HC bandpower for this sz.mat, linear units
        NMIClkHC_deltaBP_mat(:,szCount) = mean(NormMIClkHC_deltaBP,2,'omitnan'); % delta
        NMIClkHC_betaBP_mat(:,szCount) = mean(NormMIClkHC_betaBP,2,'omitnan'); % beta
        %
        % get normed Hit postictal spectrogram for this sz.mat, linear units
        NHPClkLO_SPGM_mat(:,:,szCount) = mean(NormHPClkLO_SPGM,3,'omitnan');
        NHPClkHC_SPGM_mat(:,:,szCount) = mean(NormHPClkHC_SPGM,3,'omitnan');
        % get normed Hit postictal psd for this sz.mat, linear units
        NHPClkLO_PSD_mat(:,szCount) = mean(NormHPClkLO_PSD,2,'omitnan');
        NHPClkHC_PSD_mat(:,szCount) = mean(NormHPClkHC_PSD,2,'omitnan');
        % normed Hit postictal LO bandpower for this sz.mat, linear units
        NHPClkLO_deltaBP_mat(:,szCount) = mean(NormHPClkLO_deltaBP,2,'omitnan'); % delta
        NHPClkLO_betaBP_mat(:,szCount) = mean(NormHPClkLO_betaBP,2,'omitnan'); % beta
        % normed Hit postictal HC bandpower for this sz.mat, linear units
        NHPClkHC_deltaBP_mat(:,szCount) = mean(NormHPClkHC_deltaBP,2,'omitnan'); % delta
        NHPClkHC_betaBP_mat(:,szCount) = mean(NormHPClkHC_betaBP,2,'omitnan'); % beta
        %
        % get normed Miss postictal spectrogram for this sz.mat, linear units
        NMPClkLO_SPGM_mat(:,:,szCount) = mean(NormMPClkLO_SPGM,3,'omitnan');
        NMPClkHC_SPGM_mat(:,:,szCount) = mean(NormMPClkHC_SPGM,3,'omitnan');
        % get normed Miss postictal psd for this sz.mat, linear units
        NMPClkLO_PSD_mat(:,szCount) = mean(NormMPClkLO_PSD,2,'omitnan');
        NMPClkHC_PSD_mat(:,szCount) = mean(NormMPClkHC_PSD,2,'omitnan');
        % normed Miss postictal LO bandpower for this sz.mat, linear units
        NMPClkLO_deltaBP_mat(:,szCount) = mean(NormMPClkLO_deltaBP,2,'omitnan'); % delta
        NMPClkLO_betaBP_mat(:,szCount) = mean(NormMPClkLO_betaBP,2,'omitnan'); % beta
        % normed Miss postictal HC bandpower for this sz.mat, linear units
        NMPClkHC_deltaBP_mat(:,szCount) = mean(NormMPClkHC_deltaBP,2,'omitnan'); % delta
        NMPClkHC_betaBP_mat(:,szCount) = mean(NormMPClkHC_betaBP,2,'omitnan'); % beta
        %
        % get normed Hit recovery spectrogram for this sz.mat, linear units
        NHRClkLO_SPGM_mat(:,:,szCount) = mean(NormHRClkLO_SPGM,3,'omitnan');
        NHRClkHC_SPGM_mat(:,:,szCount) = mean(NormHRClkHC_SPGM,3,'omitnan');
        % get normed Hit recovery psd for this sz.mat, linear units
        NHRClkLO_PSD_mat(:,szCount) = mean(NormHRClkLO_PSD,2,'omitnan');
        NHRClkHC_PSD_mat(:,szCount) = mean(NormHRClkHC_PSD,2,'omitnan');
        % normed Hit recovery LO bandpower for this sz.mat, linear units
        NHRClkLO_deltaBP_mat(:,szCount) = mean(NormHRClkLO_deltaBP,2,'omitnan'); % delta
        NHRClkLO_betaBP_mat(:,szCount) = mean(NormHRClkLO_betaBP,2,'omitnan'); % beta
        % normed Hit postictal HC bandpower for this sz.mat, linear units
        NHRClkHC_deltaBP_mat(:,szCount) = mean(NormHRClkHC_deltaBP,2,'omitnan'); % delta
        NHRClkHC_betaBP_mat(:,szCount) = mean(NormHRClkHC_betaBP,2,'omitnan'); % beta
        %
        % get normed Miss recovery spectrogram for this sz.mat, linear units
        NMRClkLO_SPGM_mat(:,:,szCount) = mean(NormMRClkLO_SPGM,3,'omitnan');
        NMRClkHC_SPGM_mat(:,:,szCount) = mean(NormMRClkHC_SPGM,3,'omitnan');
        % get normed Miss recovery psd for this sz.mat, linear units
        NMRClkLO_PSD_mat(:,szCount) = mean(NormMRClkLO_PSD,2,'omitnan');
        NMRClkHC_PSD_mat(:,szCount) = mean(NormMRClkHC_PSD,2,'omitnan');
        % normed Miss recovery LO bandpower for this sz.mat, linear units
        NMRClkLO_deltaBP_mat(:,szCount) = mean(NormMRClkLO_deltaBP,2,'omitnan'); % delta
        NMRClkLO_betaBP_mat(:,szCount) = mean(NormMRClkLO_betaBP,2,'omitnan'); % beta
        % normed Miss recovery HC bandpower for this sz.mat, linear units
        NMRClkHC_deltaBP_mat(:,szCount) = mean(NormMRClkHC_deltaBP,2,'omitnan'); % delta
        NMRClkHC_betaBP_mat(:,szCount) = mean(NormMRClkHC_betaBP,2,'omitnan'); % beta
    end
    %% baseline, all hit, mean across all sz mat files
    % normed baseline spectrogram, linear units
    NBClkLO_SPGM_mean = mean(NBClkLO_SPGM_mat,3,'omitnan');
    NBClkHC_SPGM_mean = mean(NBClkHC_SPGM_mat,3,'omitnan');
    % normed baseline psd, linear units
    NBClkLO_PSD_mean = mean(NBClkLO_PSD_mat,2,'omitnan');
    NBClkHC_PSD_mean = mean(NBClkHC_PSD_mat,2,'omitnan');
    % normed baseline LO bandpower, linear units
    NBClkLO_deltaBP_mean = mean(NBClkLO_deltaBP_mat,2,'omitnan'); % delta
    NBClkLO_betaBP_mean = mean(NBClkLO_betaBP_mat,2,'omitnan'); % beta
    % normed baseline HC bandpower, linear units
    NBClkHC_deltaBP_mean = mean(NBClkHC_deltaBP_mat,2,'omitnan'); % delta
    NBClkHC_betaBP_mean = mean(NBClkHC_betaBP_mat,2,'omitnan'); % beta  
    %% hit ictal, mean across all sz mat files
    % normed hit ictal spectrogram, linear units
    NHIClkLO_SPGM_mean = mean(NHIClkLO_SPGM_mat,3,'omitnan');
    NHIClkHC_SPGM_mean = mean(NHIClkHC_SPGM_mat,3,'omitnan');
    % normed hit ictal psd, linear units
    NHIClkLO_PSD_mean = mean(NHIClkLO_PSD_mat,2,'omitnan');
    NHIClkHC_PSD_mean = mean(NHIClkHC_PSD_mat,2,'omitnan');
    % normed hit ictal LO bandpower, linear units
    NHIClkLO_deltaBP_mean = mean(NHIClkLO_deltaBP_mat,2,'omitnan'); % delta
    NHIClkLO_betaBP_mean = mean(NHIClkLO_betaBP_mat,2,'omitnan'); % beta
    % normed hit ictal HC bandpower, linear units
    NHIClkHC_deltaBP_mean = mean(NHIClkHC_deltaBP_mat,2,'omitnan'); % delta
    NHIClkHC_betaBP_mean = mean(NHIClkHC_betaBP_mat,2,'omitnan'); % beta  
    %% miss ictal, mean across all sz mat files
    % normed miss ictal spectrogram, linear units
    NMIClkLO_SPGM_mean = mean(NMIClkLO_SPGM_mat,3,'omitnan');
    NMIClkHC_SPGM_mean = mean(NMIClkHC_SPGM_mat,3,'omitnan');
    % normed miss ictal psd, linear units
    NMIClkLO_PSD_mean = mean(NMIClkLO_PSD_mat,2,'omitnan');
    NMIClkHC_PSD_mean = mean(NMIClkHC_PSD_mat,2,'omitnan');
    % normed miss ictal LO bandpower, linear units
    NMIClkLO_deltaBP_mean = mean(NMIClkLO_deltaBP_mat,2,'omitnan'); % delta
    NMIClkLO_betaBP_mean = mean(NMIClkLO_betaBP_mat,2,'omitnan'); % beta
    % normed miss ictal HC bandpower, linear units
    NMIClkHC_deltaBP_mean = mean(NMIClkHC_deltaBP_mat,2,'omitnan'); % delta
    NMIClkHC_betaBP_mean = mean(NMIClkHC_betaBP_mat,2,'omitnan'); % beta
    %% Hit postictal, mean across all sz mat files
    % normed hit postictal spectrogram, linear units
    NHPClkLO_SPGM_mean = mean(NHPClkLO_SPGM_mat,3,'omitnan');
    NHPClkHC_SPGM_mean = mean(NHPClkHC_SPGM_mat,3,'omitnan');
    % normed hit postictal psd, linear units
    NHPClkLO_PSD_mean = mean(NHPClkLO_PSD_mat,2,'omitnan');
    NHPClkHC_PSD_mean = mean(NHPClkHC_PSD_mat,2,'omitnan');
    % normed hit postictal LO bandpower, linear units
    NHPClkLO_deltaBP_mean = mean(NHPClkLO_deltaBP_mat,2,'omitnan'); % delta
    NHPClkLO_betaBP_mean = mean(NHPClkLO_betaBP_mat,2,'omitnan'); % beta
    % normed hit postictal HC bandpower, linear units
    NHPClkHC_deltaBP_mean = mean(NHPClkHC_deltaBP_mat,2,'omitnan'); % delta
    NHPClkHC_betaBP_mean = mean(NHPClkHC_betaBP_mat,2,'omitnan'); % beta
    %% Miss postictal, mean across all sz mat files
    % normed miss postictal spectrogram, linear units
    NMPClkLO_SPGM_mean = mean(NMPClkLO_SPGM_mat,3,'omitnan');
    NMPClkHC_SPGM_mean = mean(NMPClkHC_SPGM_mat,3,'omitnan');
    % normed miss postictal psd, linear units
    NMPClkLO_PSD_mean = mean(NMPClkLO_PSD_mat,2,'omitnan');
    NMPClkHC_PSD_mean = mean(NMPClkHC_PSD_mat,2,'omitnan');
    % normed miss postictal LO bandpower, linear units
    NMPClkLO_deltaBP_mean = mean(NMPClkLO_deltaBP_mat,2,'omitnan'); % delta
    NMPClkLO_betaBP_mean = mean(NMPClkLO_betaBP_mat,2,'omitnan'); % beta
    % normed miss postictal HC bandpower, linear units
    NMPClkHC_deltaBP_mean = mean(NMPClkHC_deltaBP_mat,2,'omitnan'); % delta
    NMPClkHC_betaBP_mean = mean(NMPClkHC_betaBP_mat,2,'omitnan'); % beta
    %% Hit recovery, mean across all sz mat files
    % normed hit recovery spectrogram, linear units
    NHRClkLO_SPGM_mean = mean(NHRClkLO_SPGM_mat,3,'omitnan');
    NHRClkHC_SPGM_mean = mean(NHRClkHC_SPGM_mat,3,'omitnan');
    % normed hit recovery psd, linear units
    NHRClkLO_PSD_mean = mean(NHRClkLO_PSD_mat,2,'omitnan');
    NHRClkHC_PSD_mean = mean(NHRClkHC_PSD_mat,2,'omitnan');
    % normed hit recovery LO bandpower, linear units
    NHRClkLO_deltaBP_mean = mean(NHRClkLO_deltaBP_mat,2,'omitnan'); % delta
    NHRClkLO_betaBP_mean = mean(NHRClkLO_betaBP_mat,2,'omitnan'); % beta
    % normed hit recovery HC bandpower, linear units
    NHRClkHC_deltaBP_mean = mean(NHRClkHC_deltaBP_mat,2,'omitnan'); % delta
    NHRClkHC_betaBP_mean = mean(NHRClkHC_betaBP_mat,2,'omitnan'); % beta
    %% Miss recovery, mean across all sz mat files
    % normed hit recovery spectrogram, linear units
    NMRClkLO_SPGM_mean = mean(NMRClkLO_SPGM_mat,3,'omitnan');
    NMRClkHC_SPGM_mean = mean(NMRClkHC_SPGM_mat,3,'omitnan');
    % normed hit recovery psd, linear units
    NMRClkLO_PSD_mean = mean(NMRClkLO_PSD_mat,2,'omitnan');
    NMRClkHC_PSD_mean = mean(NMRClkHC_PSD_mat,2,'omitnan');
    % normed hit recovery LO bandpower, linear units
    NMRClkLO_deltaBP_mean = mean(NMRClkLO_deltaBP_mat,2,'omitnan'); % delta
    NMRClkLO_betaBP_mean = mean(NMRClkLO_betaBP_mat,2,'omitnan'); % beta
    % normed hit recovery HC bandpower, linear units
    NMRClkHC_deltaBP_mean = mean(NMRClkHC_deltaBP_mat,2,'omitnan'); % delta
    NMRClkHC_betaBP_mean = mean(NMRClkHC_betaBP_mat,2,'omitnan'); % beta
    %% output to struct
    % baseline clicks
    % spectrogram
    BClk_SPGM.BClkfy = BClkfy;
    BClk_SPGM.BClktx = BClktx;
    BClk_SPGM.NBClkLO_SPGM = NBClkLO_SPGM_mean;
    BClk_SPGM.NBClkHC_SPGM = NBClkHC_SPGM_mean;
    % PSD
    BClk_PSD.BClkfw = BClkfw;
    BClk_PSD.NBClkLO_PSD = NBClkLO_PSD_mean;
    BClk_PSD.NBClkHC_PSD = NBClkHC_PSD_mean;
    % Bandpower
    BClk_BP.NBClkLO_deltaBP = NBClkLO_deltaBP_mat; % delta band
    BClk_BP.NBClkLO_betaBP = NBClkLO_betaBP_mat; % beta band    
    BClk_BP.NBClkHC_deltaBP = NBClkHC_deltaBP_mat; % delta band
    BClk_BP.NBClkHC_betaBP = NBClkHC_betaBP_mat; % beta band
    % mean
    BClk_BP.NBClkLO_deltaBP_mean = NBClkLO_deltaBP_mean; % delta band
    BClk_BP.NBClkLO_betaBP_mean = NBClkLO_betaBP_mean; % beta band    
    BClk_BP.NBClkHC_deltaBP_mean = NBClkHC_deltaBP_mean; % delta band
    BClk_BP.NBClkHC_betaBP_mean = NBClkHC_betaBP_mean; % beta band
    %%
    % Ictal hit clicks
    % spectrogram
    IClk_SPGM.IClkfy = IClkfy;
    IClk_SPGM.IClktx = IClktx;
    IClk_SPGM.NHIClkLO_SPGM = NHIClkLO_SPGM_mean;
    IClk_SPGM.NHIClkHC_SPGM = NHIClkHC_SPGM_mean;
    % PSD
    IClk_PSD.IClkfw = IClkfw;
    IClk_PSD.NHIClkLO_PSD = NHIClkLO_PSD_mean;
    IClk_PSD.NHIClkHC_PSD = NHIClkHC_PSD_mean;
    % Bandpower
    IClk_BP.NHIClkLO_deltaBP = NHIClkLO_deltaBP_mat; % delta band
    IClk_BP.NHIClkLO_betaBP = NHIClkLO_betaBP_mat; % beta band    
    IClk_BP.NHIClkHC_deltaBP = NHIClkHC_deltaBP_mat; % delta band
    IClk_BP.NHIClkHC_betaBP = NHIClkHC_betaBP_mat; % beta band
    % mean
    IClk_BP.NHIClkLO_deltaBP_mean = NHIClkLO_deltaBP_mean; % delta band
    IClk_BP.NHIClkLO_betaBP_mean = NHIClkLO_betaBP_mean; % beta band    
    IClk_BP.NHIClkHC_deltaBP_mean = NHIClkHC_deltaBP_mean; % delta band
    IClk_BP.NHIClkHC_betaBP_mean = NHIClkHC_betaBP_mean; % beta band 
    %
    % Ictal miss clicks
    % spectrogram
    IClk_SPGM.NMIClkLO_SPGM = NMIClkLO_SPGM_mean;
    IClk_SPGM.NMIClkHC_SPGM = NMIClkHC_SPGM_mean;
    % PSD
    IClk_PSD.NMIClkLO_PSD = NMIClkLO_PSD_mean;
    IClk_PSD.NMIClkHC_PSD = NMIClkHC_PSD_mean;
    % Bandpower
    IClk_BP.NMIClkLO_deltaBP = NMIClkLO_deltaBP_mat; % delta band
    IClk_BP.NMIClkLO_betaBP = NMIClkLO_betaBP_mat; % beta band
    IClk_BP.NMIClkHC_deltaBP = NMIClkHC_deltaBP_mat; % delta band
    IClk_BP.NMIClkHC_betaBP = NMIClkHC_betaBP_mat; % beta band
    % mean
    IClk_BP.NMIClkLO_deltaBP_mean = NMIClkLO_deltaBP_mean; % delta band
    IClk_BP.NMIClkLO_betaBP_mean = NMIClkLO_betaBP_mean; % beta band
    IClk_BP.NMIClkHC_deltaBP_mean = NMIClkHC_deltaBP_mean; % delta band
    IClk_BP.NMIClkHC_betaBP_mean = NMIClkHC_betaBP_mean; % beta band
    %%
    % hit postctal clicks
    % spectrogram
    PClk_SPGM.PClkfy = PClkfy;
    PClk_SPGM.PClktx = PClktx;
    PClk_SPGM.NHPClkLO_SPGM = NHPClkLO_SPGM_mean;
    PClk_SPGM.NHPClkHC_SPGM = NHPClkHC_SPGM_mean;
    % PSD
    PClk_PSD.PClkfw = PClkfw;
    PClk_PSD.NHPClkLO_PSD = NHPClkLO_PSD_mean;
    PClk_PSD.NHPClkHC_PSD = NHPClkHC_PSD_mean;
    % Bandpower
    PClk_BP.NHPClkLO_deltaBP = NHPClkLO_deltaBP_mat; % delta band
    PClk_BP.NHPClkLO_betaBP = NHPClkLO_betaBP_mat; % beta band    
    PClk_BP.NHPClkHC_deltaBP = NHPClkHC_deltaBP_mat; % delta band
    PClk_BP.NHPClkHC_betaBP = NHPClkHC_betaBP_mat; % beta band
    % mean
    PClk_BP.NHPClkLO_deltaBP_mean = NHPClkLO_deltaBP_mean; % delta band
    PClk_BP.NHPClkLO_betaBP_mean = NHPClkLO_betaBP_mean; % beta band    
    PClk_BP.NHPClkHC_deltaBP_mean = NHPClkHC_deltaBP_mean; % delta band
    PClk_BP.NHPClkHC_betaBP_mean = NHPClkHC_betaBP_mean; % beta band
    %
    % miss postctal clicks
    % spectrogram
    PClk_SPGM.NMPClkLO_SPGM = NMPClkLO_SPGM_mean;
    PClk_SPGM.NMPClkHC_SPGM = NMPClkHC_SPGM_mean;
    % PSD
    PClk_PSD.NMPClkLO_PSD = NMPClkLO_PSD_mean;
    PClk_PSD.NMPClkHC_PSD = NMPClkHC_PSD_mean;
    % Bandpower
    PClk_BP.NMPClkLO_deltaBP = NMPClkLO_deltaBP_mat; % delta band
    PClk_BP.NMPClkLO_betaBP = NMPClkLO_betaBP_mat; % beta band    
    PClk_BP.NMPClkHC_deltaBP = NMPClkHC_deltaBP_mat; % delta band
    PClk_BP.NMPClkHC_betaBP = NMPClkHC_betaBP_mat; % beta band
    % mean
    PClk_BP.NMPClkLO_deltaBP_mean = NMPClkLO_deltaBP_mean; % delta band
    PClk_BP.NMPClkLO_betaBP_mean = NMPClkLO_betaBP_mean; % beta band    
    PClk_BP.NMPClkHC_deltaBP_mean = NMPClkHC_deltaBP_mean; % delta band
    PClk_BP.NMPClkHC_betaBP_mean = NMPClkHC_betaBP_mean; % beta band
    %%
    % hit recovery clicks
    % spectrogram
    RClk_SPGM.RClkfy = RClkfy;
    RClk_SPGM.RClktx = RClktx;
    RClk_SPGM.NHRClkLO_SPGM = NHRClkLO_SPGM_mean;
    RClk_SPGM.NHRClkHC_SPGM = NHRClkHC_SPGM_mean;
    % PSD
    RClk_PSD.RClkfw = RClkfw;
    RClk_PSD.NHRClkLO_PSD = NHRClkLO_PSD_mean;
    RClk_PSD.NHRClkHC_PSD = NHRClkHC_PSD_mean;
    % Bandpower
    RClk_BP.NHRClkLO_deltaBP = NHRClkLO_deltaBP_mat; % delta band
    RClk_BP.NHRClkLO_betaBP = NHRClkLO_betaBP_mat; % beta band    
    RClk_BP.NHRClkHC_deltaBP = NHRClkHC_deltaBP_mat; % delta band
    RClk_BP.NHRClkHC_betaBP = NHRClkHC_betaBP_mat; % beta band
    % mean
    RClk_BP.NHRClkLO_deltaBP_mean = NHRClkLO_deltaBP_mean; % delta band
    RClk_BP.NHRClkLO_betaBP_mean = NHRClkLO_betaBP_mean; % beta band    
    RClk_BP.NHRClkHC_deltaBP_mean = NHRClkHC_deltaBP_mean; % delta band
    RClk_BP.NHRClkHC_betaBP_mean = NHRClkHC_betaBP_mean; % beta band
    %
    % miss recovery clicks
    % spectrogram
    RClk_SPGM.NMRClkLO_SPGM = NMRClkLO_SPGM_mean;
    RClk_SPGM.NMRClkHC_SPGM = NMRClkHC_SPGM_mean;
    % PSD
    RClk_PSD.NMRClkLO_PSD = NMRClkLO_PSD_mean;
    RClk_PSD.NMRClkHC_PSD = NMRClkHC_PSD_mean;
    % Bandpower
    RClk_BP.NMRClkLO_deltaBP = NMRClkLO_deltaBP_mat; % delta band
    RClk_BP.NMRClkLO_betaBP = NMRClkLO_betaBP_mat; % beta band    
    RClk_BP.NMRClkHC_deltaBP = NMRClkHC_deltaBP_mat; % delta band
    RClk_BP.NMRClkHC_betaBP = NMRClkHC_betaBP_mat; % beta band
    % mean
    RClk_BP.NMRClkLO_deltaBP_mean = NMRClkLO_deltaBP_mean; % delta band
    RClk_BP.NMRClkLO_betaBP_mean = NMRClkLO_betaBP_mean; % beta band    
    RClk_BP.NMRClkHC_deltaBP_mean = NMRClkHC_deltaBP_mean; % delta band
    RClk_BP.NMRClkHC_betaBP_mean = NMRClkHC_betaBP_mean; % beta band
    %% save
    % baseline
    save(fullfile(SaveDir,'BClk_SPGM.mat'),'BClk_SPGM');
    save(fullfile(SaveDir,'BClk_PSD.mat'),'BClk_PSD');
    save(fullfile(SaveDir,'BClk_BP.mat'),'BClk_BP');
    % ictal
    save(fullfile(SaveDir,'IClk_SPGM.mat'),'IClk_SPGM');
    save(fullfile(SaveDir,'IClk_PSD.mat'),'IClk_PSD');
    save(fullfile(SaveDir,'IClk_BP.mat'),'IClk_BP');
    % postictal
    save(fullfile(SaveDir,'PClk_SPGM.mat'),'PClk_SPGM');
    save(fullfile(SaveDir,'PClk_PSD.mat'),'PClk_PSD');
    save(fullfile(SaveDir,'PClk_BP.mat'),'PClk_BP');
    % recovery
    save(fullfile(SaveDir,'RClk_SPGM.mat'),'RClk_SPGM');
    save(fullfile(SaveDir,'RClk_PSD.mat'),'RClk_PSD');
    save(fullfile(SaveDir,'RClk_BP.mat'),'RClk_BP');
    fprintf('Done!\n'); 
end