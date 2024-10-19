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
