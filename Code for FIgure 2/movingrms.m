function rms = movingrms(y, win, fs, varargin)
% function to calculate a moving rms of signal y. 
% INPUT: - y: vector or matrix to rms 
%        - win: window in seconds
%        - fs: the sampling frequency of signal 
%        - dim (optional): dimension along which to calculate rms; default
%        is 1. 
% OUTPUT: - rms: the signal after rms is calculated 
% Written by Abdo Sharaf April 18th 2021
% #########################################################################

dim = 1; 
if nargin > 3
    dim = varargin{1}; 
end

[r, ~] = size(y); 
if dim == 1 && r == 1
    y = y'; 
end

win_vals = round(win * fs);
rms = sqrt(movmean(y.^2, win_vals, dim)); 
end