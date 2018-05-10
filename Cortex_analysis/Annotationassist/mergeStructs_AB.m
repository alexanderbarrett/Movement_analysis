function [pullStruct] = mergeStructs_AB(pushStruct,pullStruct)
    

    pushNames = fieldnames(pushStruct);
    pullNames = fieldnames(pullStruct);
    % loop over all pushStruct fields, add them to pullStruct
    for name = 1:length(pushNames)
        iPushName = pushNames{name};
        % don't overwrite
        if (strmatch(iPushName,pullNames,'exact'))
            if strcmp(iPushName,'framelist_true')
                pullStruct.(iPushName) = cat(2,pullStruct.(iPushName),pushStruct.(iPushName));
            else 
                pullStruct.(iPushName) = cat(1,pullStruct.(iPushName),pushStruct.(iPushName));
            end
        else
            pullStruct.(iPushName) = pushStruct.(iPushName);
        end
    end