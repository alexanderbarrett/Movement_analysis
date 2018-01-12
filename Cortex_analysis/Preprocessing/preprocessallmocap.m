function preprocessallmocap(ratname,mocapmasterdirectory)


mocapfilestruct = loadmocapfilestruct(ratname,mocapmasterdirectory);

mocapvideodirectory = [];
%% get the desired files
descriptor_struct_1 = struct();

good_conds = setxor(fieldnames(mocapfilestruct),'mocapdir');

for cond_ind = 1:numel(good_conds)
 descriptor_struct_1.cond = good_conds{cond_ind};

for day_ind = 1:numel(mocapfilestruct.(descriptor_struct_1.cond).days);

descriptor_struct_1.day = day_ind;

for cond_ind = 1:numel( mocapfilestruct.(descriptor_struct_1.cond).day_conds{descriptor_struct_1.day})
descriptor_struct_1.tag = mocapfilestruct.(descriptor_struct_1.cond).day_conds{descriptor_struct_1.day}{cond_ind};
descriptor_struct_1.vidtag = descriptor_struct_1.tag;

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
good_inds 
if numel(good_inds)
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);

dir_base = strcat(mocapfilestruct.mocapdir,strcat(strrep(mocapfilestruct.(descriptor_struct_1.cond).days{descriptor_struct_1.day},'Generated_C3D_files\','')));
file_folders = dir(dir_base);
directories = find(cat(1,file_folders.isdir)==1);
%good_folder = [];
for mm = directories'
if (numel(strfind(file_folders(mm).name,descriptor_struct_1.vidtag)))
 %   good_folder = file_folders(mm).name;
end
end
%mocapvideodirectory =  strcat(dir_base,good_folder);




overwrite=1;
overwrite_macro = 0;
descriptor_struct_1.Vidtag_savetag = descriptor_struct_1.vidtag;
descriptor_struct_1.Condition = descriptor_struct_1.cond;
descriptor_struct_1.Days = descriptor_struct_1.day;
descriptor_struct_1.Nametag = descriptor_struct_1.vidtag;

[mocapstruct] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes,overwrite,overwrite_macro,[],[],0);
end
end
end
end
end
