% Periods based spectrogram/psd/bandpower analysis

% Edited by Jerry. 09/21/2022
function [szDur,PeriodLO_SPGM,PeriodLO_PSD,PeriodLO_BP,...
    PeriodHC_SPGM,PeriodHC_PSD,PeriodHC_BP]=...
    AnalyzePeriodFrequencyV3_220922(params,szCell,szDir,szTiming,SaveDir)
    % parse in the relevant parameters 
    sps = params.SampleRate; % 1000
    PreIctalOffset = params.PreIctalOffset*sps; % baseline 60s*1000
    IctalCutoff = params.IctalCutoff*sps; % 15s*1000 ictal cutoff
    PostIctalOffset = params.PostIctalOffset*sps; % postictal 30s*1000
    RecoveryPeriod = params.RecoveryPeriod*sps; % recovery 60s*1000
    ControlOffset = params.ControlOffset*sps; % control 15s*1000
    artifactPeriod = params.StimulationArtifactPeriod*sps;% 10s*1000
    szTimeSkip = params.SeizureOnsetTimeSkip*sps; % 2.5s*1000
    % fft number and overlap for spectrogram/psd based on period
    nff = 1000; 
    nOverlap = params.SpectrogramOverlap;  % 1
    % other parameters for spectrogram/psd
    sWin = params.SpectrogramWindow; % 1000   
    fprintf('Gathering electrophysiology information...\n'); 
    % include/exclude array for all .mat files
    toInclude = string(szCell(2:end,3));
    % freq range
    deltaf = [1 4];
    thetaf = [4 8];
    alphaf = [8 13];
    betaf = [13 30];
    %% get 60s baseline before ictal onset as reference period for psd
    szCount = 0; % number of seizures included
    for fn = 1:length(szDir)        % path to the sz .mat file 
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
        % remove 3s before AllTimesInds(2) to exclude 2s sz induction train        
        RefLO_LFP = LO_LFP.values...
            (max((AllTimesInds(1)-3000),1):AllTimesInds(2)-3000);
        RefHC_LFP = HC_LFP.values...
            (max((AllTimesInds(1)-3000),1):AllTimesInds(2)-3000);
        % compute reference spectrogram, 1000sample win
        % 1 sample overlap,1000 DFT,1000 fs,freq[Hz] shown on y.
        % output is linear unit, not dB.
        [~,REFfy,REFtx,ReferenceLO_SPGM] = spectrogram...
            (RefLO_LFP,sWin,nOverlap,nff,sps,'yaxis'); 
        [~,~,~,ReferenceHC_SPGM] = spectrogram...
            (RefHC_LFP,sWin,nOverlap,nff,sps,'yaxis');
        % actual Reference duration and Reference <= PreIctalOffset
        ReferenceDur = ...
            (AllTimesInds(2)-3000)-max((AllTimesInds(1)-3000),1)+1;
        % fill from start using nan if BaselineDur < PreIctalOffset
        if ReferenceDur < PreIctalOffset
            REFfill=nan(length(REFfy),(PreIctalOffset/sps-length(REFtx)));
            RefLO_SPGM(:,:,szCount) = cat(2,REFfill,ReferenceLO_SPGM);
            RefHC_SPGM(:,:,szCount) = cat(2,REFfill,ReferenceHC_SPGM);
        else
            RefLO_SPGM(:,:,szCount) = ...
                ReferenceLO_SPGM(:,end-PreIctalOffset/sps+1:end);
            RefHC_SPGM(:,:,szCount) = ...
                ReferenceHC_SPGM(:,end-PreIctalOffset/sps+1:end);
        end
        % compute the psd
        [RLO_PSD,REFfw] = pwelch(RefLO_LFP,sWin,nOverlap,nff,sps,'psd');
        [RHC_PSD,~] = pwelch(RefHC_LFP,sWin,nOverlap,nff,sps,'psd');
        %
        RefLO_PSD(:,szCount) = RLO_PSD;
        RefHC_PSD(:,szCount) = RHC_PSD;
        % LO bandpower, linear units
        RefLO_deltaBP(:,szCount) = bandpower(RLO_PSD,REFfw,deltaf,'psd'); 
        RefLO_thetaBP(:,szCount) = bandpower(RLO_PSD,REFfw,thetaf,'psd'); 
        RefLO_alphaBP(:,szCount) = bandpower(RLO_PSD,REFfw,alphaf,'psd'); 
        RefLO_betaBP(:,szCount) = bandpower(RLO_PSD,REFfw,betaf,'psd'); 
        % HC bandpower, linear units
        RefHC_deltaBP(:,szCount) = bandpower(RHC_PSD,REFfw,deltaf,'psd');
        RefHC_thetaBP(:,szCount) = bandpower(RHC_PSD,REFfw,thetaf,'psd');
        RefHC_alphaBP(:,szCount) = bandpower(RHC_PSD,REFfw,alphaf,'psd');
        RefHC_betaBP(:,szCount) = bandpower(RHC_PSD,REFfw,betaf,'psd');
    end
    % reference LO/HC spectrogram, linear units
    RefLO_SPGM_mean = mean(RefLO_SPGM,3,'omitnan');
    RefLO_SPGM_mean = mean(RefLO_SPGM_mean,2,'omitnan');
    RefHC_SPGM_mean = mean(RefHC_SPGM,3,'omitnan');
    RefHC_SPGM_mean = mean(RefHC_SPGM_mean,2,'omitnan');
    % reference LO/HC psd, linear units
    RefLO_PSD_mean = mean(RefLO_PSD,2,'omitnan');
    RefHC_PSD_mean = mean(RefHC_PSD,2,'omitnan');
    % reference LO bandpower, linear units
    RefLO_deltaBP_mean = mean(RefLO_deltaBP,2,'omitnan'); % delta
    RefLO_thetaBP_mean = mean(RefLO_thetaBP,2,'omitnan'); % theta
    RefLO_alphaBP_mean = mean(RefLO_alphaBP,2,'omitnan'); % alpha
    RefLO_betaBP_mean = mean(RefLO_betaBP,2,'omitnan'); % beta
    % reference HC bandpower, linear units
    RefHC_deltaBP_mean = mean(RefHC_deltaBP,2,'omitnan'); % delta
    RefHC_thetaBP_mean = mean(RefHC_thetaBP,2,'omitnan'); % theta
    RefHC_alphaBP_mean = mean(RefHC_alphaBP,2,'omitnan'); % alpha
    RefHC_betaBP_mean = mean(RefHC_betaBP,2,'omitnan'); % beta
    %% loop again for baseline,ictal,postical,recovery,and control
    szCount = 0;
    for fn = 1:length(szDir)
        % path to the sz .mat file 
        filepath = szDir{fn}; 
        disp(filepath);
        %
        if strcmp(toInclude(fn), 'Include')
            szCount = szCount + 1;  % increment the number of seizures 
        else
            continue;   % move on to next iteration 
        end
        % indeces of baseline,ictal,postical,recovery,and control
        AllTimesInds = szTiming(szCount).AllTimesInds;
        % LFP length is the sz .mat file length
        LO_LFP = szTiming(szCount).LO_LFP; 
        HC_LFP = szTiming(szCount).HC_LFP; 
        hasControl = szTiming(szCount).hasControl;
        szDur{szCount} = AllTimesInds(4)-AllTimesInds(3)+1;
        %% remove artifact following sz/control stimulation using polyfit
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
        %% LFP spectrogram/psd for control,baseline,ictal,postical,recovery
        % baseline LFP spectrogram/psd
        % read baseline LFP, remove 3s before AllTimesInds(2)
        BaselineLO_LFP = LO_LFP.values...
            (max((AllTimesInds(1)-3000),1):(AllTimesInds(2)-3000));
        BaselineHC_LFP = HC_LFP.values...
            (max((AllTimesInds(1)-3000),1):(AllTimesInds(2)-3000));
        % compute baseline spectrogram, 1000sample win
        % 1 sample overlap,1000 DFT,1000 fs,freq[Hz] shown on y.
        [~,bfy,btx,BaselineLO_SPGM] = spectrogram...
            (BaselineLO_LFP,sWin,nOverlap,nff,sps,'yaxis'); 
        [~,~,~,BaselineHC_SPGM] = spectrogram...
            (BaselineHC_LFP,sWin,nOverlap,nff,sps,'yaxis');
        % normalization using reference
        NBaselineLO_SPGM = BaselineLO_SPGM./RefLO_SPGM_mean;
        NBaselineHC_SPGM = BaselineHC_SPGM./RefHC_SPGM_mean;
        % actual baseline duration and BaselineDur <= PreIctalOffset
        BaselineDur = ...
            (AllTimesInds(2)-3000)-max((AllTimesInds(1)-3000),1)+1;
        % fill from start using nan if BaselineDur < PreIctalOffset
        if BaselineDur < PreIctalOffset
            bfill = nan(length(bfy),(PreIctalOffset/sps-length(btx)));
            BLO_SPGM(:,:,szCount) = cat(2,bfill,BaselineLO_SPGM);
            BHC_SPGM(:,:,szCount) = cat(2,bfill,BaselineHC_SPGM);
            NormBLO_SPGM(:,:,szCount) = cat(2,bfill,NBaselineLO_SPGM);
            NormBHC_SPGM(:,:,szCount) = cat(2,bfill,NBaselineHC_SPGM);
        else
            BLO_SPGM(:,:,szCount) = ...
                BaselineLO_SPGM(:,end-PreIctalOffset/sps+1:end);
            BHC_SPGM(:,:,szCount) = ...
                BaselineHC_SPGM(:,end-PreIctalOffset/sps+1:end);
            NormBLO_SPGM(:,:,szCount) = ...
                NBaselineLO_SPGM(:,end-PreIctalOffset/sps+1:end);
            NormBHC_SPGM(:,:,szCount) = ...
                NBaselineHC_SPGM(:,end-PreIctalOffset/sps+1:end);
        end
        % compute the psd
        [BaselineLO_PSD,bfw] = pwelch...
            (BaselineLO_LFP,sWin,nOverlap,nff,sps,'psd');
        [BaselineHC_PSD,~] = pwelch...
            (BaselineHC_LFP,sWin,nOverlap,nff,sps,'psd');
        %
        BLO_PSD(:,szCount) = BaselineLO_PSD;
        BHC_PSD(:,szCount) = BaselineHC_PSD;
        % normalization using reference
        NormBLO_PSD(:,szCount) = BaselineLO_PSD./RefLO_PSD_mean;
        NormBHC_PSD(:,szCount) = BaselineHC_PSD./RefHC_PSD_mean;
        % bandpower LO
        BLO_deltaBP(:,szCount) = ...
            bandpower(BaselineLO_PSD,bfw,deltaf,'psd'); % delta
        BLO_thetaBP(:,szCount) = ...
            bandpower(BaselineLO_PSD,bfw,thetaf,'psd'); % theta
        BLO_alphaBP(:,szCount) = ...
            bandpower(BaselineLO_PSD,bfw,alphaf,'psd'); % alpha
        BLO_betaBP(:,szCount) = ...
            bandpower(BaselineLO_PSD,bfw,betaf,'psd'); % beta
        % normalization using reference
        NormBLO_deltaBP(:,szCount) = ...
            bandpower(BaselineLO_PSD,bfw,deltaf,'psd')/RefLO_deltaBP_mean;
        NormBLO_thetaBP(:,szCount) = ...
            bandpower(BaselineLO_PSD,bfw,thetaf,'psd')/RefLO_thetaBP_mean;
        NormBLO_alphaBP(:,szCount) = ...
            bandpower(BaselineLO_PSD,bfw,alphaf,'psd')/RefLO_alphaBP_mean;
        NormBLO_betaBP(:,szCount) = ...
            bandpower(BaselineLO_PSD,bfw,betaf,'psd')/RefLO_betaBP_mean;
        % bandpower HC
        BHC_deltaBP(:,szCount) = ...
            bandpower(BaselineHC_PSD,bfw,deltaf,'psd'); % delta
        BHC_thetaBP(:,szCount) = ...
            bandpower(BaselineHC_PSD,bfw,thetaf,'psd'); % theta
        BHC_alphaBP(:,szCount) = ...
            bandpower(BaselineHC_PSD,bfw,alphaf,'psd'); % alpha
        BHC_betaBP(:,szCount) = ...
            bandpower(BaselineHC_PSD,bfw,betaf,'psd'); % beta
        % normalization using reference
        NormBHC_deltaBP(:,szCount) = ...
            bandpower(BaselineHC_PSD,bfw,deltaf,'psd')/RefHC_deltaBP_mean;
        NormBHC_thetaBP(:,szCount) = ...
            bandpower(BaselineHC_PSD,bfw,thetaf,'psd')/RefHC_thetaBP_mean;
        NormBHC_alphaBP(:,szCount) = ...
            bandpower(BaselineHC_PSD,bfw,alphaf,'psd')/RefHC_alphaBP_mean;
        NormBHC_betaBP(:,szCount) = ...
            bandpower(BaselineHC_PSD,bfw,betaf,'psd')/RefHC_betaBP_mean;
        %% ictal LFP spectrogram/psd, remember to remove szTimeSkip part
        IctalStart = AllTimesInds(3)+szTimeSkip; % ictal onset
        % read ictal LFP
        IctalLO_LFP = LO_LFP.values(IctalStart:AllTimesInds(4));
        IctalHC_LFP = HC_LFP.values(IctalStart:AllTimesInds(4));
        % remove short sz
        TempDur = AllTimesInds(4)-IctalStart+1;
        if TempDur<sWin
            % spectrogram
            ILO_SPGM(:,:,szCount) = nan(length(bfy),IctalCutoff/sps);
            IHC_SPGM(:,:,szCount) = nan(length(bfy),IctalCutoff/sps);
            NormILO_SPGM(:,:,szCount) = nan(length(bfy),IctalCutoff/sps);
            NormIHC_SPGM(:,:,szCount) = nan(length(bfy),IctalCutoff/sps);
            % PSD
            ILO_PSD(:,szCount) = nan(length(bfy),1);
            IHC_PSD(:,szCount) = nan(length(bfy),1);
            NormILO_PSD(:,szCount) = nan(length(bfy),1);
            NormIHC_PSD(:,szCount) = nan(length(bfy),1);
            % LO Bandpower
            ILO_deltaBP(:,szCount) = nan;
            ILO_thetaBP(:,szCount) = nan;
            ILO_alphaBP(:,szCount) = nan;
            ILO_betaBP(:,szCount) = nan;
            NormILO_deltaBP(:,szCount) = nan;
            NormILO_thetaBP(:,szCount) = nan;
            NormILO_alphaBP(:,szCount) = nan;
            NormILO_betaBP(:,szCount) = nan;
            % HC Bandpower
            IHC_deltaBP(:,szCount) = nan;
            IHC_thetaBP(:,szCount) = nan;
            IHC_alphaBP(:,szCount) = nan;
            IHC_betaBP(:,szCount) = nan;
            NormIHC_deltaBP(:,szCount) = nan;
            NormIHC_thetaBP(:,szCount) = nan;
            NormIHC_alphaBP(:,szCount) = nan;
            NormIHC_betaBP(:,szCount) = nan;
        else
            % compute ictal spectrogram, 1000sample win
            % 1 sample overlap,1000 DFT,1000 fs,freq[Hz] shown on y.
            [~,ify,itx,IctalLO_SPGM] = spectrogram...
                (IctalLO_LFP,sWin,nOverlap,nff,sps,'yaxis'); 
            [~,~,~,IctalHC_SPGM] = spectrogram...
                (IctalHC_LFP,sWin,nOverlap,nff,sps,'yaxis');
            % normalization using reference
            NIctalLO_SPGM = IctalLO_SPGM./RefLO_SPGM_mean;
            NIctalHC_SPGM = IctalHC_SPGM./RefHC_SPGM_mean;
            % actual ictal duration
            IctalDur = AllTimesInds(4)-IctalStart+1; 
            % fill with nan to the end if ictal duration less than IctalCutoff
            if IctalDur < IctalCutoff
                ifill = nan(length(ify),(IctalCutoff/sps-length(itx)));
                ILO_SPGM(:,:,szCount) = cat(2,IctalLO_SPGM,ifill);
                IHC_SPGM(:,:,szCount) = cat(2,IctalHC_SPGM,ifill);
                NormILO_SPGM(:,:,szCount) = cat(2,NIctalLO_SPGM,ifill);
                NormIHC_SPGM(:,:,szCount) = cat(2,NIctalHC_SPGM,ifill);
            % get the first IctalCutoff/sps
            else
                ILO_SPGM(:,:,szCount) = IctalLO_SPGM(:,1:IctalCutoff/sps);
                IHC_SPGM(:,:,szCount) = IctalHC_SPGM(:,1:IctalCutoff/sps);
                NormILO_SPGM(:,:,szCount) = ...
                    NIctalLO_SPGM(:,1:IctalCutoff/sps);
                NormIHC_SPGM(:,:,szCount) = ...
                    NIctalHC_SPGM(:,1:IctalCutoff/sps);
            end        
            % compute the psd
            [IctalLO_PSD,ifw] = pwelch...
                (IctalLO_LFP,sWin,nOverlap,nff,sps,'psd');
            [IctalHC_PSD,~] = pwelch...
                (IctalHC_LFP,sWin,nOverlap,nff,sps,'psd');
            %
            ILO_PSD(:,szCount) = IctalLO_PSD;
            IHC_PSD(:,szCount) = IctalHC_PSD;
            % normalization using reference
            NormILO_PSD(:,szCount) = IctalLO_PSD./RefLO_PSD_mean;
            NormIHC_PSD(:,szCount) = IctalHC_PSD./RefHC_PSD_mean;
            % bandpower LO
            ILO_deltaBP(:,szCount) = ...
                bandpower(IctalLO_PSD,ifw,deltaf,'psd'); % delta
            ILO_thetaBP(:,szCount) = ...
                bandpower(IctalLO_PSD,ifw,thetaf,'psd'); % theta
            ILO_alphaBP(:,szCount) = ...
                bandpower(IctalLO_PSD,ifw,alphaf,'psd'); % alpha
            ILO_betaBP(:,szCount) = ...
                bandpower(IctalLO_PSD,ifw,betaf,'psd'); % beta
            % normalization using reference
            NormILO_deltaBP(:,szCount) = ...
                bandpower(IctalLO_PSD,ifw,deltaf,'psd')/RefLO_deltaBP_mean;
            NormILO_thetaBP(:,szCount) = ...
                bandpower(IctalLO_PSD,ifw,thetaf,'psd')/RefLO_thetaBP_mean;
            NormILO_alphaBP(:,szCount) = ...
                bandpower(IctalLO_PSD,ifw,alphaf,'psd')/RefLO_alphaBP_mean;
            NormILO_betaBP(:,szCount) = ...
                bandpower(IctalLO_PSD,ifw,betaf,'psd')/RefLO_betaBP_mean;
            % bandpower HC
            IHC_deltaBP(:,szCount) = ...
                bandpower(IctalHC_PSD,ifw,deltaf,'psd'); % delta
            IHC_thetaBP(:,szCount) = ...
                bandpower(IctalHC_PSD,ifw,thetaf,'psd'); % theta
            IHC_alphaBP(:,szCount) = ...
                bandpower(IctalHC_PSD,ifw,alphaf,'psd'); % alpha
            IHC_betaBP(:,szCount) = ...
                bandpower(IctalHC_PSD,ifw,betaf,'psd'); % beta
            % normalization using reference
            NormIHC_deltaBP(:,szCount) = ...
                bandpower(IctalHC_PSD,ifw,deltaf,'psd')/RefHC_deltaBP_mean;
            NormIHC_thetaBP(:,szCount) = ...
                bandpower(IctalHC_PSD,ifw,thetaf,'psd')/RefHC_thetaBP_mean;
            NormIHC_alphaBP(:,szCount) = ...
                bandpower(IctalHC_PSD,ifw,alphaf,'psd')/RefHC_alphaBP_mean;
            NormIHC_betaBP(:,szCount) = ...
                bandpower(IctalHC_PSD,ifw,betaf,'psd')/RefHC_betaBP_mean;
        end  
        %% postictal LFP spectrogram/psd
        % read postictal LFP
        PostictalLO_LFP = LO_LFP.values(AllTimesInds(5):AllTimesInds(6));
        PostictalHC_LFP = HC_LFP.values(AllTimesInds(5):AllTimesInds(6));
        % compute postictal spectrogram, 1000sample win
        % 1 sample overlap,1000 DFT,1000 fs,freq[Hz] shown on y.
        [~,pfy,ptx,PostictalLO_SPGM] = spectrogram...
            (PostictalLO_LFP,sWin,nOverlap,nff,sps,'yaxis'); 
        [~,~,~,PostictalHC_SPGM] = spectrogram...
            (PostictalHC_LFP,sWin,nOverlap,nff,sps,'yaxis');
        % normalization using reference
        NPostictalLO_SPGM = PostictalLO_SPGM./RefLO_SPGM_mean;
        NPostictalHC_SPGM = PostictalHC_SPGM./RefHC_SPGM_mean;
        % actual postictal duration
        PostictalDur = AllTimesInds(6)-AllTimesInds(5)+1; 
        % fill with nan if postictal duration less than PostIctalOffset
        if PostictalDur < PostIctalOffset
            pfill = nan(length(pfy),(PostIctalOffset/sps-length(ptx)));
            PLO_SPGM(:,:,szCount) = cat(2,PostictalLO_SPGM,pfill);
            PHC_SPGM(:,:,szCount) = cat(2,PostictalHC_SPGM,pfill);
            NormPLO_SPGM(:,:,szCount) = cat(2,NPostictalLO_SPGM,pfill);
            NormPHC_SPGM(:,:,szCount) = cat(2,NPostictalHC_SPGM,pfill);
        % get the first PostIctalOffset/sps
        else
            PLO_SPGM(:,:,szCount) = ...
                PostictalLO_SPGM(:,1:PostIctalOffset/sps);
            PHC_SPGM(:,:,szCount) = ...
                PostictalHC_SPGM(:,1:PostIctalOffset/sps);
            NormPLO_SPGM(:,:,szCount) = ...
                NPostictalLO_SPGM(:,1:PostIctalOffset/sps);
            NormPHC_SPGM(:,:,szCount) = ...
                NPostictalHC_SPGM(:,1:PostIctalOffset/sps);
        end
        % compute the psd
        [PostictalLO_PSD,pfw] = pwelch...
            (PostictalLO_LFP,sWin,nOverlap,nff,sps,'psd');
        [PostictalHC_PSD,~] = pwelch...
            (PostictalHC_LFP,sWin,nOverlap,nff,sps,'psd');
        %
        PLO_PSD(:,szCount) = PostictalLO_PSD;
        PHC_PSD(:,szCount) = PostictalHC_PSD;
        % normalization using reference
        NormPLO_PSD(:,szCount) = PostictalLO_PSD./RefLO_PSD_mean;
        NormPHC_PSD(:,szCount) = PostictalHC_PSD./RefHC_PSD_mean;
        % bandpower LO
        PLO_deltaBP(:,szCount) = ...
            bandpower(PostictalLO_PSD,pfw,deltaf,'psd'); % delta
        PLO_thetaBP(:,szCount) = ...
            bandpower(PostictalLO_PSD,pfw,thetaf,'psd'); % theta
        PLO_alphaBP(:,szCount) = ...
            bandpower(PostictalLO_PSD,pfw,alphaf,'psd'); % alpha
        PLO_betaBP(:,szCount) = ...
            bandpower(PostictalLO_PSD,pfw,betaf,'psd'); % beta
        % normalization using reference
        NormPLO_deltaBP(:,szCount) = ...
            bandpower(PostictalLO_PSD,pfw,deltaf,'psd')/RefLO_deltaBP_mean;
        NormPLO_thetaBP(:,szCount) = ...
            bandpower(PostictalLO_PSD,pfw,thetaf,'psd')/RefLO_thetaBP_mean;
        NormPLO_alphaBP(:,szCount) = ...
            bandpower(PostictalLO_PSD,pfw,alphaf,'psd')/RefLO_alphaBP_mean;
        NormPLO_betaBP(:,szCount) = ...
            bandpower(PostictalLO_PSD,pfw,betaf,'psd')/RefLO_betaBP_mean;
        % bandpower HC
        PHC_deltaBP(:,szCount) = ...
            bandpower(PostictalHC_PSD,pfw,deltaf,'psd'); % delta
        PHC_thetaBP(:,szCount) = ...
            bandpower(PostictalHC_PSD,pfw,thetaf,'psd'); % theta
        PHC_alphaBP(:,szCount) = ...
            bandpower(PostictalHC_PSD,pfw,alphaf,'psd'); % alpha
        PHC_betaBP(:,szCount) = ...
            bandpower(PostictalHC_PSD,pfw,betaf,'psd'); % beta
        % normalization using reference
        NormPHC_deltaBP(:,szCount) = ...
            bandpower(PostictalHC_PSD,pfw,deltaf,'psd')/RefHC_deltaBP_mean;
        NormPHC_thetaBP(:,szCount) = ...
            bandpower(PostictalHC_PSD,pfw,thetaf,'psd')/RefHC_thetaBP_mean;
        NormPHC_alphaBP(:,szCount) = ...
            bandpower(PostictalHC_PSD,pfw,alphaf,'psd')/RefHC_alphaBP_mean;
        NormPHC_betaBP(:,szCount) = ...
            bandpower(PostictalHC_PSD,pfw,betaf,'psd')/RefHC_betaBP_mean;
        %% recovery LFP spectrogram/psd, 
        % 1:length(LO_LFP.values) is the index of sz.mat file
        if length(LO_LFP.values)<=AllTimesInds(7) % No recovery, fill nan
            RecoveryLO_SPGM = nan(length(bfy),RecoveryPeriod/sps);
            RecoveryHC_SPGM = nan(length(bfy),RecoveryPeriod/sps);
            NormRecoveryLO_SPGM = nan(length(bfy),RecoveryPeriod/sps);
            NormRecoveryHC_SPGM = nan(length(bfy),RecoveryPeriod/sps);
        % Part of the recovery, fill the remainings with nans
        elseif length(LO_LFP.values)<=AllTimesInds(8)
            % actual recovery duration
            RecoveryDur = length(LO_LFP.values)-AllTimesInds(7)+1;
            % read in recovery LFP
            RecoveryLO_LFP = LO_LFP.values...
                (AllTimesInds(7):AllTimesInds(7)+RecoveryDur-1);
            RecoveryHC_LFP = HC_LFP.values...
                (AllTimesInds(7):AllTimesInds(7)+RecoveryDur-1);
            % compute recovery spectrogram, 1000sample win
            % 1 sample overlap,1000 DFT,1000 fs,freq[Hz] shown on y.
            [~,rfy,rtx,RecoveryLO_SPGM] = spectrogram...
                (RecoveryLO_LFP,sWin,nOverlap,nff,sps,'yaxis'); 
            [~,~,~,RecoveryHC_SPGM] = spectrogram...
                (RecoveryHC_LFP,sWin,nOverlap,nff,sps,'yaxis');
            % normalization using reference
            NRecoveryLO_SPGM = RecoveryLO_SPGM./RefLO_SPGM_mean;
            NRecoveryHC_SPGM = RecoveryHC_SPGM./RefHC_SPGM_mean;
            % fill with nan if RecoveryDur less than RecoveryPeriod
            if RecoveryDur < RecoveryPeriod
                rfill = nan(length(rfy),(RecoveryPeriod/sps-length(rtx)));
                RecoveryLO_SPGM = cat(2,RecoveryLO_SPGM,rfill);
                RecoveryHC_SPGM = cat(2,RecoveryHC_SPGM,rfill);
                NRecoveryLO_SPGM = cat(2,NRecoveryLO_SPGM,rfill);
                NRecoveryHC_SPGM = cat(2,NRecoveryHC_SPGM,rfill);
            % get the first RecoveryPeriod/sps
            else
                RecoveryLO_SPGM = RecoveryLO_SPGM(:,1:RecoveryPeriod/sps);
                RecoveryHC_SPGM = RecoveryHC_SPGM(:,1:RecoveryPeriod/sps);
                NRecoveryLO_SPGM=NRecoveryLO_SPGM(:,1:RecoveryPeriod/sps);
                NRecoveryHC_SPGM=NRecoveryHC_SPGM(:,1:RecoveryPeriod/sps);
            end
        else
            % actual recovery duration
            RecoveryDur = AllTimesInds(8)-AllTimesInds(7)+1;
            % read in recovery LFP
            RecoveryLO_LFP=LO_LFP.values(AllTimesInds(7):AllTimesInds(8));
            RecoveryHC_LFP=HC_LFP.values(AllTimesInds(7):AllTimesInds(8));
            % compute recovery spectrogram, 1000sample win
            % 1 sample overlap,1000 DFT,1000 fs,freq[Hz] shown on y.
            [~,rfy,rtx,RecoveryLO_SPGM] = spectrogram...
                (RecoveryLO_LFP,sWin,nOverlap,nff,sps,'yaxis'); 
            [~,~,~,RecoveryHC_SPGM] = spectrogram...
                (RecoveryHC_LFP,sWin,nOverlap,nff,sps,'yaxis');
            % normalization using reference
            NRecoveryLO_SPGM = RecoveryLO_SPGM./RefLO_SPGM_mean;
            NRecoveryHC_SPGM = RecoveryHC_SPGM./RefHC_SPGM_mean;
            % fill with nan if RecoveryDur less than RecoveryPeriod
            if RecoveryDur < RecoveryPeriod
                rfill=nan(length(rfy),(RecoveryPeriod/sps-length(rtx)));
                RecoveryLO_SPGM = cat(2,RecoveryLO_SPGM,rfill);
                RecoveryHC_SPGM = cat(2,RecoveryHC_SPGM,rfill);
                NRecoveryLO_SPGM = cat(2,NRecoveryLO_SPGM,rfill);
                NRecoveryHC_SPGM = cat(2,NRecoveryHC_SPGM,rfill);
            % get the first RecoveryPeriod/sps
            else
                RecoveryLO_SPGM = RecoveryLO_SPGM(:,1:RecoveryPeriod/sps);
                RecoveryHC_SPGM = RecoveryHC_SPGM(:,1:RecoveryPeriod/sps);
                NRecoveryLO_SPGM=NRecoveryLO_SPGM(:,1:RecoveryPeriod/sps);
                NRecoveryHC_SPGM=NRecoveryHC_SPGM(:,1:RecoveryPeriod/sps);
            end
        end
        RLO_SPGM(:,:,szCount) = RecoveryLO_SPGM;
        RHC_SPGM(:,:,szCount) = RecoveryHC_SPGM;
        NormRLO_SPGM(:,:,szCount) = NRecoveryLO_SPGM;
        NormRHC_SPGM(:,:,szCount) = NRecoveryHC_SPGM;
        % compute the psd
        [RecoveryLO_PSD,rfw] = pwelch...
            (RecoveryLO_LFP,sWin,nOverlap,nff,sps,'psd');
        [RecoveryHC_PSD,~] = pwelch...
            (RecoveryHC_LFP,sWin,nOverlap,nff,sps,'psd');
        RLO_PSD(:,szCount) = RecoveryLO_PSD;
        RHC_PSD(:,szCount) = RecoveryHC_PSD;
        % normalization using reference
        NormRLO_PSD(:,szCount) = RecoveryLO_PSD./RefLO_PSD_mean;
        NormRHC_PSD(:,szCount) = RecoveryHC_PSD./RefHC_PSD_mean;
        % bandpower LO
        RLO_deltaBP(:,szCount) = ...
            bandpower(RecoveryLO_PSD,rfw,deltaf,'psd'); % delta
        RLO_thetaBP(:,szCount) = ...
            bandpower(RecoveryLO_PSD,rfw,deltaf,'psd'); % theta
        RLO_alphaBP(:,szCount) = ...
            bandpower(RecoveryLO_PSD,rfw,alphaf,'psd'); % alpha
        RLO_betaBP(:,szCount) = ...
            bandpower(RecoveryLO_PSD,rfw,betaf,'psd'); % beta
        % normalization using reference
        NormRLO_deltaBP(:,szCount) = ...
            bandpower(RecoveryLO_PSD,rfw,deltaf,'psd')/RefLO_deltaBP_mean;
        NormRLO_thetaBP(:,szCount) = ...
            bandpower(RecoveryLO_PSD,rfw,thetaf,'psd')/RefLO_thetaBP_mean;
        NormRLO_alphaBP(:,szCount) = ...
            bandpower(RecoveryLO_PSD,rfw,alphaf,'psd')/RefLO_alphaBP_mean;
        NormRLO_betaBP(:,szCount) = ...
            bandpower(RecoveryLO_PSD,rfw,betaf,'psd')/RefLO_betaBP_mean;
        % bandpower HC
        RHC_deltaBP(:,szCount) = ...
            bandpower(RecoveryHC_PSD,rfw,deltaf,'psd'); % delta
        RHC_thetaBP(:,szCount) = ...
            bandpower(RecoveryHC_PSD,rfw,thetaf,'psd'); % theta
        RHC_alphaBP(:,szCount) = ...
            bandpower(RecoveryHC_PSD,rfw,alphaf,'psd'); % alpha
        RHC_betaBP(:,szCount) = ...
            bandpower(RecoveryHC_PSD,rfw,betaf,'psd'); % beta
        % normalization using reference
        NormRHC_deltaBP(:,szCount) = ...
            bandpower(RecoveryHC_PSD,rfw,deltaf,'psd')/RefHC_deltaBP_mean;
        NormRHC_thetaBP(:,szCount) = ...
            bandpower(RecoveryHC_PSD,rfw,thetaf,'psd')/RefHC_thetaBP_mean;
        NormRHC_alphaBP(:,szCount) = ...
            bandpower(RecoveryHC_PSD,rfw,alphaf,'psd')/RefHC_alphaBP_mean;
        NormRHC_betaBP(:,szCount) = ...
            bandpower(RecoveryHC_PSD,rfw,betaf,'psd')/RefHC_betaBP_mean;
        %% control LFP spectrogram/psd, 
        if hasControl % has AllTimesInds(9),(10)
            % read control LFP
            ControlLO_LFP=LO_LFP.values(AllTimesInds(9):AllTimesInds(10));
            ControlHC_LFP=HC_LFP.values(AllTimesInds(9):AllTimesInds(10));
            % remove short control
            TempDur = AllTimesInds(10)-AllTimesInds(9)+1;
            if TempDur<sWin
                % spectrogram
                CLO_SPGM(:,:,szCount) = nan(length(bfy),ControlOffset/sps);
                CHC_SPGM(:,:,szCount) = nan(length(bfy),ControlOffset/sps);
                NormCLO_SPGM(:,:,szCount) = ...
                    nan(length(bfy),ControlOffset/sps);
                NormCHC_SPGM(:,:,szCount) = ...
                    nan(length(bfy),ControlOffset/sps);
                % PSD
                CLO_PSD(:,szCount) = nan(length(bfy),1);
                CHC_PSD(:,szCount) = nan(length(bfy),1);
                NormCLO_PSD(:,szCount) = nan(length(bfy),1);
                NormCHC_PSD(:,szCount) = nan(length(bfy),1);
                % LO Bandpower
                CLO_deltaBP(:,szCount) = nan;
                CLO_thetaBP(:,szCount) = nan;
                CLO_alphaBP(:,szCount) = nan;
                CLO_betaBP(:,szCount) = nan;
                NormCLO_deltaBP(:,szCount) = nan;
                NormCLO_thetaBP(:,szCount) = nan;
                NormCLO_alphaBP(:,szCount) = nan;
                NormCLO_betaBP(:,szCount) = nan;
                % HC Bandpower
                CHC_deltaBP(:,szCount) = nan;
                CHC_thetaBP(:,szCount) = nan;
                CHC_alphaBP(:,szCount) = nan;
                CHC_betaBP(:,szCount) = nan;
                NormCHC_deltaBP(:,szCount) = nan;
                NormCHC_thetaBP(:,szCount) = nan;
                NormCHC_alphaBP(:,szCount) = nan;
                NormCHC_betaBP(:,szCount) = nan;
            else
                % compute control spectrogram (psd), 1000sample win
                % 1 sample overlap,1000 DFT,1000 fs,freq[Hz] shown on y.
                [~,cfy,ctx,ControlLO_SPGM] = spectrogram...
                    (ControlLO_LFP,sWin,nOverlap,nff,sps,'yaxis'); 
                [~,~,~,ControlHC_SPGM] = spectrogram...
                    (ControlHC_LFP,sWin,nOverlap,nff,sps,'yaxis');
                % normalization using reference
                NControlLO_SPGM = ControlLO_SPGM./RefLO_SPGM_mean;
                NControlHC_SPGM = ControlHC_SPGM./RefHC_SPGM_mean;
                % actual control duration
                ControlDur = AllTimesInds(10)-AllTimesInds(9)+1;
                % fill the end using nan
                if ControlDur < ControlOffset
                    cfill=nan(length(cfy),(ControlOffset/sps-length(ctx)));
                    CLO_SPGM(:,:,szCount) = cat(2,ControlLO_SPGM,cfill);
                    CHC_SPGM(:,:,szCount) = cat(2,ControlHC_SPGM,cfill);
                    NormCLO_SPGM(:,:,szCount)=cat(2,NControlLO_SPGM,cfill);
                    NormCHC_SPGM(:,:,szCount)=cat(2,NControlHC_SPGM,cfill);
                % get the first ControlOff/sps
                else
                    CLO_SPGM(:,:,szCount) = ...
                        ControlLO_SPGM(:,1:ControlOffset/sps);
                    CHC_SPGM(:,:,szCount) = ...
                        ControlHC_SPGM(:,1:ControlOffset/sps);
                    NormCLO_SPGM(:,:,szCount) = ...
                        NControlLO_SPGM(:,1:ControlOffset/sps);
                    NormCHC_SPGM(:,:,szCount) = ...
                        NControlHC_SPGM(:,1:ControlOffset/sps);
                end
                % compute the psd
                [ControlLO_PSD,cfw] = pwelch...
                    (ControlLO_LFP,sWin,nOverlap,nff,sps,'psd');
                [ControlHC_PSD,~] = pwelch...
                    (ControlHC_LFP,sWin,nOverlap,nff,sps,'psd');
                CLO_PSD(:,szCount) = ControlLO_PSD;
                CHC_PSD(:,szCount) = ControlHC_PSD;
                % normalization using reference
                NormCLO_PSD(:,szCount) = ControlLO_PSD./RefLO_PSD_mean;
                NormCHC_PSD(:,szCount) = ControlHC_PSD./RefHC_PSD_mean;
                % bandpower LO
                CLO_deltaBP(:,szCount) = ...
                    bandpower(ControlLO_PSD,cfw,deltaf,'psd'); % delta
                CLO_thetaBP(:,szCount) = ...
                    bandpower(ControlLO_PSD,cfw,thetaf,'psd'); % theta
                CLO_alphaBP(:,szCount) = ...
                    bandpower(ControlLO_PSD,cfw,alphaf,'psd'); % alpha
                CLO_betaBP(:,szCount) = ...
                    bandpower(ControlLO_PSD,cfw,betaf,'psd'); % beta
                % normalization using reference
                NormCLO_deltaBP(:,szCount) = bandpower...
                    (ControlLO_PSD,cfw,deltaf,'psd')/RefLO_deltaBP_mean;
                NormCLO_thetaBP(:,szCount) = bandpower...
                    (ControlLO_PSD,cfw,thetaf,'psd')/RefLO_thetaBP_mean;
                NormCLO_alphaBP(:,szCount) = bandpower...
                    (ControlLO_PSD,cfw,alphaf,'psd')/RefLO_alphaBP_mean;
                NormCLO_betaBP(:,szCount) = bandpower...
                    (ControlLO_PSD,cfw,betaf,'psd')/RefLO_betaBP_mean;
                % bandpower HC
                CHC_deltaBP(:,szCount) = ...
                    bandpower(ControlHC_PSD,cfw,deltaf,'psd'); % delta
                CHC_thetaBP(:,szCount) = ...
                    bandpower(ControlHC_PSD,cfw,thetaf,'psd'); % theta
                CHC_alphaBP(:,szCount) = ...
                    bandpower(ControlHC_PSD,cfw,alphaf,'psd'); % alpha
                CHC_betaBP(:,szCount) = ...
                    bandpower(ControlHC_PSD,cfw,betaf,'psd'); % beta
                % normalization using reference
                NormCHC_deltaBP(:,szCount) = bandpower...
                    (ControlHC_PSD,cfw,deltaf,'psd')/RefHC_deltaBP_mean;
                NormCHC_thetaBP(:,szCount) = bandpower...
                    (ControlHC_PSD,cfw,thetaf,'psd')/RefHC_thetaBP_mean;
                NormCHC_alphaBP(:,szCount) = bandpower...
                    (ControlHC_PSD,cfw,alphaf,'psd')/RefHC_alphaBP_mean;
                NormCHC_betaBP(:,szCount) = bandpower...
                    (ControlHC_PSD,cfw,betaf,'psd')/RefHC_betaBP_mean;
            end
        end
    end
    %% Reference LO
    PeriodLO_SPGM.REFfy = REFfy;
    PeriodLO_SPGM.REFtx = REFtx;
    PeriodLO_SPGM.ReferenceLO_SPGM = RefLO_SPGM;
    PeriodLO_PSD.bfwelch = REFfw;
    PeriodLO_PSD.ReferenceLO_PSD = RefLO_PSD;
    PeriodLO_BP.ReferenceDeltaBP = RefLO_deltaBP; % delta band
    PeriodLO_BP.ReferenceThetaBP = RefLO_thetaBP; % theta band
    PeriodLO_BP.ReferenceAlphaBP = RefLO_alphaBP; % alpha band
    PeriodLO_BP.ReferenceBetaBP = RefLO_betaBP; % beta band
    % Baseline HC
    PeriodHC_SPGM.bfy = REFfy;
    PeriodHC_SPGM.btx = REFtx;
    PeriodHC_SPGM.ReferenceHC_SPGM = RefHC_SPGM;
    PeriodHC_PSD.bfwelch = REFfw;
    PeriodHC_PSD.ReferenceHC_PSD = RefHC_PSD;
    PeriodHC_BP.ReferenceDeltaBP = RefHC_deltaBP; % delta band
    PeriodHC_BP.ReferenceThetaBP = RefHC_thetaBP; % theta band
    PeriodHC_BP.ReferenceAlphaBP = RefHC_alphaBP; % alpha band
    PeriodHC_BP.ReferenceBetaBP = RefHC_betaBP; % beta band
    %% Baseline LO
    PeriodLO_SPGM.bfy = bfy;
    PeriodLO_SPGM.btx = btx;
    PeriodLO_SPGM.BaselineLO_SPGM = BLO_SPGM;
    PeriodLO_SPGM.NormBaselineLO_SPGM = NormBLO_SPGM;
    PeriodLO_PSD.bfwelch = bfw;
    PeriodLO_PSD.BaselineLO_PSD = BLO_PSD;
    PeriodLO_PSD.NormBaselineLO_PSD = NormBLO_PSD;
    PeriodLO_BP.BaselineDeltaBP = BLO_deltaBP; % delta band
    PeriodLO_BP.BaselineThetaBP = BLO_thetaBP; % theta band
    PeriodLO_BP.BaselineAlphaBP = BLO_alphaBP; % alpha band
    PeriodLO_BP.BaselineBetaBP = BLO_betaBP; % beta band
    PeriodLO_BP.NormBaselineDeltaBP = NormBLO_deltaBP; % Norm delta band
    PeriodLO_BP.NormBaselineThetaBP = NormBLO_thetaBP; % Norm theta band
    PeriodLO_BP.NormBaselineAlphaBP = NormBLO_alphaBP; % Norm alpha band
    PeriodLO_BP.NormBaselineBetaBP = NormBLO_betaBP; % Norm beta band
    % Baseline HC
    PeriodHC_SPGM.bfy = bfy;
    PeriodHC_SPGM.btx = btx;
    PeriodHC_SPGM.BaselineHC_SPGM = BHC_SPGM;
    PeriodHC_SPGM.NormBaselineHC_SPGM = NormBHC_SPGM;
    PeriodHC_PSD.bfwelch = bfw;
    PeriodHC_PSD.BaselineHC_PSD = BHC_PSD;
    PeriodHC_PSD.NormBaselineHC_PSD = NormBHC_PSD;
    PeriodHC_BP.BaselineDeltaBP = BHC_deltaBP; % delta band
    PeriodHC_BP.BaselineThetaBP = BHC_thetaBP; % theta band
    PeriodHC_BP.BaselineAlphaBP = BHC_alphaBP; % alpha band
    PeriodHC_BP.BaselineBetaBP = BHC_betaBP; % beta band
    PeriodHC_BP.NormBaselineDeltaBP = NormBHC_deltaBP; % Norm delta band
    PeriodHC_BP.NormBaselineThetaBP = NormBHC_thetaBP; % Norm theta band
    PeriodHC_BP.NormBaselineAlphaBP = NormBHC_alphaBP; % Norm alpha band
    PeriodHC_BP.NormBaselineBetaBP = NormBHC_betaBP; % Norm beta band
    %% Ictal LO
    PeriodLO_SPGM.ify = ify;
    PeriodLO_SPGM.itx = itx;
    PeriodLO_SPGM.IctalLO_SPGM = ILO_SPGM;
    PeriodLO_SPGM.NormIctalLO_SPGM = NormILO_SPGM;
    PeriodLO_PSD.ifwelch = ifw;
    PeriodLO_PSD.IctalLO_PSD = ILO_PSD;
    PeriodLO_PSD.NormIctalLO_PSD = NormILO_PSD;
    PeriodLO_BP.IctalDeltaBP = ILO_deltaBP; % delta band
    PeriodLO_BP.IctalThetaBP = ILO_thetaBP; % theta band
    PeriodLO_BP.IctalAlphaBP = ILO_alphaBP; % alpha band
    PeriodLO_BP.IctalBetaBP = ILO_betaBP; % beta band
    PeriodLO_BP.NormIctalDeltaBP = NormILO_deltaBP; % Norm delta band
    PeriodLO_BP.NormIctalThetaBP = NormILO_thetaBP; % Norm theta band
    PeriodLO_BP.NormIctalAlphaBP = NormILO_alphaBP; % Norm alpha band
    PeriodLO_BP.NormIctalBetaBP = NormILO_betaBP; % Norm beta band
    % Ictal HC
    PeriodHC_SPGM.ify = ify;
    PeriodHC_SPGM.itx = itx;
    PeriodHC_SPGM.IctalHC_SPGM = IHC_SPGM;
    PeriodHC_SPGM.NormIctalHC_SPGM = NormIHC_SPGM;
    PeriodHC_PSD.ifwelch = ifw;
    PeriodHC_PSD.IctalHC_PSD = IHC_PSD;
    PeriodHC_PSD.NormIctalHC_PSD = NormIHC_PSD;
    PeriodHC_BP.IctalDeltaBP = IHC_deltaBP; % delta band
    PeriodHC_BP.IctalThetaBP = IHC_thetaBP; % theta band
    PeriodHC_BP.IctalAlphaBP = IHC_alphaBP; % alpha band
    PeriodHC_BP.IctalBetaBP = IHC_betaBP; % beta band
    PeriodHC_BP.NormIctalDeltaBP = NormIHC_deltaBP; % Norm delta band
    PeriodHC_BP.NormIctalThetaBP = NormIHC_thetaBP; % Norm theta band
    PeriodHC_BP.NormIctalAlphaBP = NormIHC_alphaBP; % Norm alpha band
    PeriodHC_BP.NormIctalBetaBP = NormIHC_betaBP; % Norm beta band
    %% Postictal LO
    PeriodLO_SPGM.pfy = pfy;
    PeriodLO_SPGM.ptx = ptx;
    PeriodLO_SPGM.PostictalLO_SPGM = PLO_SPGM;
    PeriodLO_SPGM.NormPostictalLO_SPGM = NormPLO_SPGM;
    PeriodLO_PSD.pfwelch = pfw;
    PeriodLO_PSD.PostictalLO_PSD = PLO_PSD;
    PeriodLO_PSD.NormPostictalLO_PSD = NormPLO_PSD;
    PeriodLO_BP.PostictalDeltaBP = PLO_deltaBP; % delta band
    PeriodLO_BP.PostictalThetaBP = PLO_thetaBP; % theta band
    PeriodLO_BP.PostictalAlphaBP = PLO_alphaBP; % alpha band
    PeriodLO_BP.PostictalBetaBP = PLO_betaBP; % beta band
    PeriodLO_BP.NormPostictalDeltaBP = NormPLO_deltaBP; % Norm delta band
    PeriodLO_BP.NormPostictalThetaBP = NormPLO_thetaBP; % Norm theta band
    PeriodLO_BP.NormPostictalAlphaBP = NormPLO_alphaBP; % Norm alpha band
    PeriodLO_BP.NormPostictalBetaBP = NormPLO_betaBP; % Norm beta band
    % Postictal HC
    PeriodHC_SPGM.pfy = pfy;
    PeriodHC_SPGM.ptx = ptx;
    PeriodHC_SPGM.PostictalHC_SPGM = PHC_SPGM;
    PeriodHC_SPGM.NormPostictalHC_SPGM = NormPHC_SPGM;
    PeriodHC_PSD.pfwelch = pfw;
    PeriodHC_PSD.PostictalHC_PSD = PHC_PSD;
    PeriodHC_PSD.NormPostictalHC_PSD = NormPHC_PSD;
    PeriodHC_BP.PostictalDeltaBP = PHC_deltaBP; % delta band
    PeriodHC_BP.PostictalThetaBP = PHC_thetaBP; % theta band
    PeriodHC_BP.PostictalAlphaBP = PHC_alphaBP; % alpha band
    PeriodHC_BP.PostictalBetaBP = PHC_betaBP; % beta band
    PeriodHC_BP.NormPostictalDeltaBP = NormPHC_deltaBP; % Norm delta band
    PeriodHC_BP.NormPostictalThetaBP = NormPHC_thetaBP; % Norm theta band
    PeriodHC_BP.NormPostictalAlphaBP = NormPHC_alphaBP; % Norm alpha band
    PeriodHC_BP.NormPostictalBetaBP = NormPHC_betaBP; % Norm beta band
    %% Recovery LO
    PeriodLO_SPGM.rfy = rfy;
    PeriodLO_SPGM.rtx = rtx;
    PeriodLO_SPGM.RecoveryLO_SPGM = RLO_SPGM;
    PeriodLO_SPGM.NormRecoveryLO_SPGM = NormRLO_SPGM;
    PeriodLO_PSD.rfwelch = rfw;
    PeriodLO_PSD.RecoveryLO_PSD = RLO_PSD;
    PeriodLO_PSD.NormRecoveryLO_PSD = NormRLO_PSD;
    PeriodLO_BP.RecoveryDeltaBP = RLO_deltaBP; % delta band
    PeriodLO_BP.RecoveryThetaBP = RLO_thetaBP; % theta band
    PeriodLO_BP.RecoveryAlphaBP = RLO_alphaBP; % alpha band
    PeriodLO_BP.RecoveryBetaBP = RLO_betaBP; % beta band
    PeriodLO_BP.NormRecoveryDeltaBP = NormRLO_deltaBP; % Norm delta band
    PeriodLO_BP.NormRecoveryThetaBP = NormRLO_thetaBP; % Norm theta band
    PeriodLO_BP.NormRecoveryAlphaBP = NormRLO_alphaBP; % Norm alpha band
    PeriodLO_BP.NormRecoveryBetaBP = NormRLO_betaBP; % Norm beta band
    % Postictal HC
    PeriodHC_SPGM.rfy = rfy;
    PeriodHC_SPGM.rtx = rtx;
    PeriodHC_SPGM.RecoveryHC_SPGM = RHC_SPGM;
    PeriodHC_SPGM.NormRecoveryHC_SPGM = NormRHC_SPGM;
    PeriodHC_PSD.rfwelch = rfw;
    PeriodHC_PSD.RecoveryHC_PSD = RHC_PSD;
    PeriodHC_PSD.NormRecoveryHC_PSD = NormRHC_PSD;
    PeriodHC_BP.RecoveryDeltaBP = RHC_deltaBP; % delta band
    PeriodHC_BP.RecoveryThetaBP = RHC_thetaBP; % theta band
    PeriodHC_BP.RecoveryAlphaBP = RHC_alphaBP; % alpha band
    PeriodHC_BP.RecoveryBetaBP = RHC_betaBP; % beta band
    PeriodHC_BP.NormRecoveryDeltaBP = NormRHC_deltaBP; % Norm delta band
    PeriodHC_BP.NormRecoveryThetaBP = NormRHC_thetaBP; % Norm theta band
    PeriodHC_BP.NormRecoveryAlphaBP = NormRHC_alphaBP; % Norm alpha band
    PeriodHC_BP.NormRecoveryBetaBP = NormRHC_betaBP; % Norm beta band
    %% control LO
    PeriodLO_SPGM.cfy = cfy;
    PeriodLO_SPGM.ctx = ctx;
    PeriodLO_SPGM.ControlLO_SPGM = CLO_SPGM; % spectrogram
    PeriodLO_SPGM.NormControlLO_SPGM = NormCLO_SPGM; % norm spectrogram
    PeriodLO_PSD.cfwelch = cfw;
    PeriodLO_PSD.ControlLO_PSD = CLO_PSD; % psd
    PeriodLO_PSD.NormControlLO_PSD = NormCLO_PSD; % Norm psd
    PeriodLO_BP.ControlDeltaBP = CLO_deltaBP; % delta band
    PeriodLO_BP.ControlThetaBP = CLO_thetaBP; % theta band
    PeriodLO_BP.ControlAlphaBP = CLO_alphaBP; % alpha band
    PeriodLO_BP.ControlBetaBP = CLO_betaBP; % beta band
    PeriodLO_BP.NormControlDeltaBP = NormCLO_deltaBP; % Norm delta band
    PeriodLO_BP.NormControlThetaBP = NormCLO_thetaBP; % Norm theta band
    PeriodLO_BP.NormControlAlphaBP = NormCLO_alphaBP; % Norm alpha band
    PeriodLO_BP.NormControlBetaBP = NormCLO_betaBP; % Norm beta band
    % control HC
    PeriodHC_SPGM.cfy = cfy;
    PeriodHC_SPGM.ctx = ctx;
    PeriodHC_SPGM.ControlHC_SPGM = CHC_SPGM;
    PeriodLO_SPGM.NormControlHC_SPGM = NormCHC_SPGM; % norm spectrogram
    PeriodHC_PSD.cfwelch = cfw;
    PeriodHC_PSD.ControlHC_PSD = CHC_PSD;
    PeriodHC_PSD.NormControlHC_PSD = NormCHC_PSD; % Norm psd
    PeriodHC_BP.ControlDeltaBP = CHC_deltaBP; % delta band
    PeriodHC_BP.ControlThetaBP = CHC_thetaBP; % theta band
    PeriodHC_BP.ControlAlphaBP = CHC_alphaBP; % alpha band
    PeriodHC_BP.ControlBetaBP = CHC_betaBP; % beta band
    PeriodHC_BP.NormControlDeltaBP = NormCHC_deltaBP; % Norm delta band
    PeriodHC_BP.NormControlThetaBP = NormCHC_thetaBP; % Norm theta band
    PeriodHC_BP.NormControlAlphaBP = NormCHC_alphaBP; % Norm alpha band
    PeriodHC_BP.NormControlBetaBP = NormCHC_betaBP; % Norm beta band
    %% save
    save(fullfile(SaveDir,'PeriodLO_SPGM.mat'),'PeriodLO_SPGM');
    save(fullfile(SaveDir,'PeriodLO_PSD.mat'),'PeriodLO_PSD');
    save(fullfile(SaveDir,'PeriodLO_BP.mat'),'PeriodLO_BP');
    save(fullfile(SaveDir,'PeriodHC_SPGM.mat'),'PeriodHC_SPGM');
    save(fullfile(SaveDir,'PeriodHC_PSD.mat'),'PeriodHC_PSD');
    save(fullfile(SaveDir,'PeriodHC_BP.mat'),'PeriodHC_BP');
    fprintf('Done!\n'); 
end