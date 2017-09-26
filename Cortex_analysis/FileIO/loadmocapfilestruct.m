
function mocapfilestruct = loadmocapfilestruct(ratname,mocapmasterdirectory)

load(strcat(mocapmasterdirectory,ratname,filesep,'mocapfilestruct_',ratname,'_.mat'))

filestruct_conds = fieldnames(mocapfilestruct);
for ll = 1:numel(filestruct_conds)
     if isfield(mocapfilestruct.(filestruct_conds{ll}),'days')
         fprintf('For type %s nhours %f nframes %f \n',(filestruct_conds{ll}),...
             sum(cellfun(@sum,mocapfilestruct.(filestruct_conds{ll}).numframes))./(300*3600),...
             sum(cellfun(@sum,mocapfilestruct.(filestruct_conds{ll}).numframes)));
     end
end

end
%mocapfilestruct
%save(strcat(mocapfilestruct.mocapdir,'mocapfilestruct_',ratname,'_.mat'),'mocapfilestruct');


