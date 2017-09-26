function [foldernums,nameFolds] = getfoldernumbers(videopathFolder)

     
     d = dir(videopathFolder);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
     foldernums = cellfun(@str2num,nameFolds);
end