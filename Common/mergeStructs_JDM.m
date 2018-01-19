function [pullStruct] = mergeStructs_JDM(pushStruct,pullStruct)
    % copies fields in pushStruct into pullStruct, if there is an overlap in field names, pushStruct overwrites pullStruct unless specified otherwise.
    %
    % [pullStruct] = mergeStructs(pushStruct,pullStruct,overwritePullFields)
    %
    % biafra ahanonu
    % started: 2014.02.12
    %
    % inputs
    % 	pushStruct -
    % 	pullStruct -
    % 	overwritePullFields - 1 = overwrite pullStruct fields with pushStruct, 0 = don't overwrite.
    % outputs
    % 	pullStruct - structure with pushStructs added

    % changelog
    %
    % TODO
    %

    pushNames = fieldnames(pushStruct);
    pullNames = fieldnames(pullStruct);
    % loop over all pushStruct fields, add them to pushStruct
    for name = 1:length(pushNames)
        iPushName = pushNames{name};
        % don't overwrite
    	if (strmatch(iPushName,pullNames,'exact'))
            pullStruct.(iPushName) = cat(2,pullStruct.(iPushName),pushStruct.(iPushName));
        else
      pullStruct.(iPushName) = pushStruct.(iPushName);
        end
    end