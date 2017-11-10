function [descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files(rat,save_tag,mocapmasterdirectory)


switch rat
    case 'Vicon8'
mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);
    case 'JDM25'
mocapfilestruct = loadmocapfilestruct('JDM25',mocapmasterdirectory);
end
mocapvideodirectory = [];
%% get the desired files
descriptor_struct_1 = struct();

switch save_tag
    case 'Vicon8_prelesion'
   descriptor_struct_1.day = 5;
descriptor_struct_1.tag = 'overnight1';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'overnight';

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);

dir_base = strcat(mocapfilestruct.mocapdir,strcat(strrep(mocapfilestruct.PreLesion.days{5},'Generated_C3D_files\','')));
file_folders = dir(dir_base);
directories = find(cat(1,file_folders.isdir)==1);
good_folder = [];
for mm = directories'
if (numel(strfind(file_folders(mm).name,descriptor_struct_1.vidtag)))
    good_folder = file_folders(mm).name;
end
end
mocapvideodirectory =  strcat(dir_base,good_folder);




case 'Vicon8_videotest'
   descriptor_struct_1.day = 5;
descriptor_struct_1.tag = 'social';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'social';

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);


dir_base = strcat(mocapfilestruct.mocapdir,strcat(strrep(mocapfilestruct.PreLesion.days{5},'Generated_C3D_files\','')));
file_folders = dir(dir_base);
directories = find(cat(1,file_folders.isdir)==1);
good_folder = [];
for mm = directories'
if (numel(strfind(file_folders(mm).name,descriptor_struct_1.vidtag)))
    good_folder = file_folders(mm).name;
end
end
mocapvideodirectory =  strcat(dir_base,good_folder);


    case 'Vicon8_caff'
descriptor_struct_1.tag = 'caff';
descriptor_struct_1.cond = 'PreLesion';
   descriptor_struct_1.day = 6;
descriptor_struct_1.vidtag = 'caff';

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
 mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);

 

dir_base = strcat(mocapfilestruct.mocapdir,strcat(strrep(mocapfilestruct.(descriptor_struct_1.cond).days{descriptor_struct_1.day},'Generated_C3D_files\','')));
file_folders = dir(dir_base);
directories = find(cat(1,file_folders.isdir)==1);
good_folder = [];
for mm = directories'
if (numel(strfind(file_folders(mm).name,descriptor_struct_1.vidtag)))
    good_folder = file_folders(mm).name;
end
end
mocapvideodirectory =  strcat(dir_base,good_folder);
 
 
    case 'Vicon8_amph'
        descriptor_struct_1.tag = 'amph';
descriptor_struct_1.cond = 'PreLesion';
   descriptor_struct_1.day = 7;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds); 
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);
 descriptor_struct_1.vidtag = 'amph';

 case 'Vicon8_dlslesion_early'
       descriptor_struct_1.tag = 'overnight_two1';
descriptor_struct_1.cond = 'UniLesion';
   descriptor_struct_1.day = 2;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);
descriptor_struct_1.vidtag = 'dlsearly';


 case 'Vicon8_prelesion_long'
     descriptor_struct_1.cond = 'PreLesion';
 descriptor_struct_1.day = 5;
 descriptor_struct_1.vidtag = 'long';
descriptor_struct_1.tag = 'long';

mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [5:7]
  descriptor_struct_1.day = day_find;

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(:));
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(:));

 end
     
 
  case 'Vicon8_lesion_long'
     descriptor_struct_1.cond = 'UniLesion';
 descriptor_struct_1.day = 5;
 descriptor_struct_1.vidtag = 'postlong';
descriptor_struct_1.tag = 'postlong';

mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [3:7]
  descriptor_struct_1.day = day_find;

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(:));
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(:));

 end
     
 case 'Vicon8_dlslesion_late'
    descriptor_struct_1.tag = 'overnight';
descriptor_struct_1.cond = 'UniLesion';
   descriptor_struct_1.day = 7;
   descriptor_struct_1.vidtag = 'post_late';

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds(1:10));
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);



      case 'Vicon8_veh'
              descriptor_struct_1.tag = 'veh';
descriptor_struct_1.cond = 'PreLesion';
   descriptor_struct_1.day = 5;
   descriptor_struct_1.vidtag = 'vehicle';

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);
 
           case 'JDM25pre'
descriptor_struct_1.tag = 'overnight1';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM25_pre';

   descriptor_struct_1.day = 6;
   good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);

 case 'JDM25post'
descriptor_struct_1.tag = 'overnight1';
descriptor_struct_1.cond = 'UniLesion';
descriptor_struct_1.vidtag = 'JDM25_post';

   descriptor_struct_1.day = 6;
   good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);


end

