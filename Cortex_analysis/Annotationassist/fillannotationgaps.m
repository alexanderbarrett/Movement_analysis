function [fulloutput,indivbouts] = fillannotationgaps(inputarray,numframes)

inputarray = unique(sort(inputarray,'ASCEND'));

diff_input = diff(inputarray);
fulloutput = [];
%define borders
diff_input(diff_input>numframes) = 0;
diff_input = [numframes diff_input];
outputstruct = bwconncomp(diff_input);
indivbouts = cell(1,numel(outputstruct.PixelIdxList));

for mm = 1:numel(outputstruct.PixelIdxList)
   fulloutput = cat(2,fulloutput,inputarray( outputstruct.PixelIdxList{mm}(1)):inputarray(outputstruct.PixelIdxList{mm}(end)));
   indivbouts{mm} = inputarray( outputstruct.PixelIdxList{mm}(1)):inputarray(outputstruct.PixelIdxList{mm}(end));
end


end