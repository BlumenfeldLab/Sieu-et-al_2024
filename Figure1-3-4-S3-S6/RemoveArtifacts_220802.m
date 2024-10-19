% Edited by Jerry. 08/02/2022
function processedData = RemoveArtifacts_220802(rawData,polyOrder)
    % *********************************************************************
    % this function process a signal to remove polynomial artifacts. 
    % Inputs: 
    % (1) rawData: raw data with artifact 
    % (2) polyOrder: order of the polynomial for polyfit
    % Outputs: 
    % (1) processedData: processed data without artifact 
    % *********************************************************************
    % make it a column vector to avoid inconsistent dimensions errors 
    if isrow(rawData)
        rawData = rawData';
    end
    % fit a polynomial to the original signal 
    d = 1:length(rawData);
    % mu(1) is mean(rawData), and mu(2) is std(rawData).
    [p,~,mu] = polyfit(d',rawData,polyOrder);
    % evaluate the fitted polynomial
    f_y = polyval(p,d,[],mu);
    % process
    processedData = rawData-f_y';
end