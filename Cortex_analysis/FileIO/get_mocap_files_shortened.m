function [descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory)


switch descriptor_struct.Rat
    case 'Vicon3'
mocapfilestruct = loadmocapfilestruct('Vicon3',mocapmasterdirectory);
%mocapfilestruct.mocapdir = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\';

    case 'Vicon8'
mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);
    case 'JDM25'
mocapfilestruct = loadmocapfilestruct('JDM25',mocapmasterdirectory);
   case 'JDM32'
mocapfilestruct = loadmocapfilestruct('JDM32',mocapmasterdirectory);
  case 'JDM21'
mocapfilestruct = loadmocapfilestruct('JDM21',mocapmasterdirectory);
end
mocapvideodirectory = [];
%% get the desired files


mocapfilearray=[];
 mocapfiletimes = [];
 for day_find = descriptor_struct.Days
descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct.Condition).mocapfiles{day_find},descriptor_struct.Searchtag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct.Condition).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct.Condition).threshcrossings{day_find}(good_inds);
switch descriptor_struct.Task
    case 0
nontaskfiles = find(mocapthreshcrossings<descriptor_struct.TaskThreshold);
case 1
nontaskfiles = find(mocapthreshcrossings>descriptor_struct.TaskThreshold);
case 2
nontaskfiles = good_inds;
end
good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct.Condition).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct.Condition).mocapdatecreate{day_find}(good_inds)');
 end
 
  goodinds = union(find(datetime(mocapfiletimes).Hour>descriptor_struct.HourRange_early) , find(datetime(mocapfiletimes).Hour<descriptor_struct.HourRange_late));
 numsample = min(numel(goodinds),descriptor_struct.NumFiles);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
 

  %% do the video part here
%   
% %videofiles
% dir_base = strcat(mocapfilestruct.mocapdir,strcat(strrep(mocapfilestruct.PreLesion.days{5},'Generated_C3D_files\','')));
% file_folders = dir(dir_base);
% directories = find(cat(1,file_folders.isdir)==1);
% good_folder = [];
% for mm = directories'
% if (numel(strfind(file_folders(mm).name,descriptor_struct_1.vidtag)))
%     good_folder = file_folders(mm).name;
% end
% end
% mocapvideodirectory =  strcat(dir_base,good_folder);
  %------------------------



end