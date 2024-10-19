% finish spectrogram,psd,and bandpower calcalation
% Edited by Jerry. 10/11/2022
function [Clkfy,Clktx,Clk_SPGM,Clkfw,Clk_PSD,Clk_deltaBP,Clk_betaBP]=...
    GetFrequencyAnalysis_221011(params,ClkLFP)
    % get parameters
    % spectrogram and fft parameters
    sps = params.SampleRate; % sample rate 
    specNfft = params.SpectrogramNfft; 
    specWindow = params.SpectrogramWindow;
    overlap = params.EventSpectrogramOverlap; 
    % freq range
    deltaf = [1 4];
    betaf = [13 30];
    % spectrogram, linear units
    [~,Clkfy,Clktx,Clk_SPGM] = ...
        spectrogram(ClkLFP,specWindow,overlap,specNfft,sps,'yaxis');
    % psd, linear units
    [Clk_PSD,Clkfw] = pwelch(ClkLFP,specWindow,overlap,specNfft,sps,'psd');
    % bandpower, linear units
    % delta
    Clk_deltaBP = bandpower(Clk_PSD,Clkfw,deltaf,'psd');
    % beta
    Clk_betaBP = bandpower(Clk_PSD,Clkfw,betaf,'psd');
end