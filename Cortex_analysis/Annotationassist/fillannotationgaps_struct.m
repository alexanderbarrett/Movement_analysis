function [fulloutput_sctuct,indivbouts_struct] = fillannotationgaps_struct(inputarray_struct,numframes)
fulloutput_sctuct = inputarray_struct;
indivbouts_struct = inputarray_struct;

fieldnamesloop = fieldnames(inputarray_struct);
for mm = 1:numel(fieldnamesloop)
    if (numel(inputarray_struct.(fieldnamesloop{mm})>5))
        [fulloutput_sctuct.(fieldnamesloop{mm}),indivbouts_struct.(fieldnamesloop{mm})] =...
            fillannotationgaps(inputarray_struct.(fieldnamesloop{mm}),numframes);
    end



end
