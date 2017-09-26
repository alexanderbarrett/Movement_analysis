function savemocappreprocessed()

mocapfilestruct = loadmocapfilestruct('Vicon8');

for dayprocess = 5:7
%get the desired files
goodinds = find(cellfun(@numel,strfind(mocapfilestruct.PreLesion.mocapfiles{dayprocess},'.c3d')));
mocapfilearray = mocapfilestruct.PreLesion.mocapfiles{dayprocess}(goodinds);

preproc_save_directory = strcat(mocapfilestruct.mocapdir,strrep(mocapfilestruct.PreLesion.days{dayprocess},'Generated_C3D_files\','Preprocessed\'));
filebase_directory = strcat(mocapfilestruct.mocapdir,mocapfilestruct.PreLesion.days{dayprocess});
mkdir(preproc_save_directory);
for mm =1:numel(mocapfilearray)
%sort them
filearrayhere{1} = mocapfilearray{mm};

[mocap_struct] =...
    preprocess_mocap_data(filearrayhere{1},mocapfilestruct.mocapdir);
fps = 300;

filename = strcat(preproc_save_directory,strrep(strrep(filearrayhere{1},filebase_directory,''),'.c3d','.mat'));
save(filename,'mocap_struct','-v7.3')

end
end
end