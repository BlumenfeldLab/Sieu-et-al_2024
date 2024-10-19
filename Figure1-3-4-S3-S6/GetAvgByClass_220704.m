% Edited by Jerry. 07/04/2022
function avgs = GetAvgByClass_220704(A,classIDs)
    %********************************************************************
    % Average an array by 'type'.  
    % Class could be seizure: 
    % Array elements corresponding to certain seizures
    % and each element corresponding to one seizure will be averaged.
    % Same thing if the class was animal.
    % Inputs:
    % (1) A: array to be averaged
    % (2) classIDs: identifiers for each element in A which will be 
    %     used to group elements in the same class together
    % Outputs:
    % (1) avgs: averaged array across elements of the same classes
    %********************************************************************
    % check if the inputs are the same size
    if size(A, 1) ~= length(classIDs)
        error('A and classIDs must be the same size');
    end
    % initialize ouputs
    avgs = [];
    % get uniqueIDs for the classes
    % if the classIDs is a string array, get numeric IDs for each class.
    % A class's ID is the index 
    % of the first occurence of that class in the classIDs array
    if isstring(classIDs)
        [~, unqIDs, ~] = unique(classIDs);
    else
        unqIDs = unique(classIDs);
    end
    for i = 1:length(unqIDs)
        avgs(i, 1) = unqIDs(i);   % ID corresponding to this average
        % number of elements in this class
        if isstring(classIDs)
            avgs(i,2) = sum(classIDs==classIDs(unqIDs(i)));
            % average
            avgs(i,3:2+size(A,2)) = ...
                mean(A(classIDs==classIDs(unqIDs(i)),:),'omitnan'); 
        else
            avgs(i,2) = sum(classIDs==unqIDs(i));
            % average
            avgs(i,3:2+size(A,2)) = ...
                mean(A(classIDs==unqIDs(i), :),'omitnan'); 
        end
    end
end
