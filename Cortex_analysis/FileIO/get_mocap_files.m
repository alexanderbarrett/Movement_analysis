function [descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files(rat,save_tag,mocapmasterdirectory)
%% INPUTS
% rat:
%save_tag: 
%mocapmasterdir
%% OUTPUTS

switch rat
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
descriptor_struct_1 = struct();


standard_num_files = 24;
standard_range_late = 17;
standard_range_early = 7;
drug_numfiles = 6;
std_prelesion_days = 3:6;


switch save_tag
    case 'Vicon8_prelesion'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'overnight';

  mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [2:5]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );






%% get video folders
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


case 'Vicon8_prelesion2'
   descriptor_struct_1.day = 5;
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'overnight';

  mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [5]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );

case 'Vicon8_prelesion_task'
   descriptor_struct_1.day = 5;
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'task';

  mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [2:6]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings>10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );




case 'Vicon8_videotest'
   descriptor_struct_1.day = 5;
descriptor_struct_1.tag = 'social';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'social';

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);


%videofiles
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

%-----------------------------
    case 'Vicon8_caff'
descriptor_struct_1.tag = 'caff';
descriptor_struct_1.cond = 'PreLesion';
   descriptor_struct_1.day = 6;
descriptor_struct_1.vidtag = 'caff';

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds(1:drug_numfiles));
 mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds(1:drug_numfiles));

 
%% videofiles
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
 


    case 'Vicon3_htr'
descriptor_struct_1.tag = 'amph';
descriptor_struct_1.cond = 'eighteenmarker';
   descriptor_struct_1.day = 2;
descriptor_struct_1.vidtag = 'amph';

%good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));

mocapfilearray = {'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\20170721\Generated_C3D_files\Vicon3_recording_amphetamine_vid1_nolj.c3d',...
    'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\20170721\Generated_C3D_files\Vicon3_recording_amphetamine_vid2_nolj.c3d',...
    'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\20170721\Generated_C3D_files\Vicon3_recording_amphetamine_vid3_nolj.c3d',...
    'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\20170721\Generated_C3D_files\Vicon3_recording_amphetamine_vid4_nolj.c3d',...
    'Y:\Jesse\Data\Motionanalysis_captures\Vicon3\20170721\Generated_C3D_files\Vicon3_recording_amphetamine_vid5_nolj.c3d',...
    };



good_inds = 1:5;
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
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds(1:drug_numfiles)); 
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds(1:drug_numfiles));
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
descriptor_struct_1.tag = 'Rec';

mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [5:7]
  descriptor_struct_1.day = day_find;

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(:));
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(:));

 end
     
 %% --------------------------------------------------------------
  case 'Vicon8_lesion_long'
     descriptor_struct_1.cond = 'UniLesion';
 descriptor_struct_1.day = 5;
 descriptor_struct_1.vidtag = 'postlong';
descriptor_struct_1.tag = 'Rec';

mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [6:8]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );


     
  %% --------------------------------------------------------------
 case 'Vicon8_dlslesion_late'
    descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'UniLesion';
   descriptor_struct_1.day = 7;
   descriptor_struct_1.vidtag = 'post_late';

  mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [7:8]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );

  
  
  
  %%--------------------------------------------
   case 'Vicon8_dlslesion_late'
    descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'UniLesion';
   descriptor_struct_1.day = 7;
   descriptor_struct_1.vidtag = 'post_late';

  mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [7:8]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );

  
  %% --------------------------------------------------------------
 case 'Vicon8_dlslesion_task'
    descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'UniLesion';
   descriptor_struct_1.day = 7;
   descriptor_struct_1.vidtag = 'post_late2_task';

  mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [6:8]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings>10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
  


 %% --------------------------------------------------------------
      case 'Vicon8_veh'
              descriptor_struct_1.tag = 'veh';
descriptor_struct_1.cond = 'PreLesion';
   descriptor_struct_1.day = 5;
   descriptor_struct_1.vidtag = 'vehicle';

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds);
 

 %% --------------------------------------------------------------
           case 'JDM25pre'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM25_pre';

  %JDM25 pre has a wonky 10x analog issue
taskthreshold = 10000;

  
  mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [6:7]
     if (strcmp(descriptor_struct_1.cond ,'PreLesion')&& day_find > 6)
         taskthreshold = 100000;
     end
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<taskthreshold);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );



 %% --------------------------------------------------------------
  case 'JDM25pre_long'
descriptor_struct_1.tag = 'overnight';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM25_pre_long';

   descriptor_struct_1.day = 6:8;
   mocapfilearray=[];
 mocapfiletimes = [];
 
   
 for day_find = [6:8]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);


mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 %%---------------------------------------------
 case 'JDM25pre2'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM25_pre2';

   descriptor_struct_1.day = 4:5;
   mocapfilearray=[];
 mocapfiletimes = [];
 
   
 for day_find = [4:5]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);


mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
 %% --------------------------------------------------------------
           case 'JDM25pre_task'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM25_pre_task';

  %JDM25 pre has a wonky 10x analog issue
taskthreshold = 10000;

  
  mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [4:7]
     if (strcmp(descriptor_struct_1.cond ,'PreLesion')&& day_find > 6)
         taskthreshold = 100000;
     end
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings>taskthreshold);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );

%% --------------------------------------------------------------
           case 'JDM25pre_all'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM25_pre_all';

  %JDM25 pre has a wonky 10x analog issue
taskthreshold = 10000;

  
  mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [4:7]
     if (strcmp(descriptor_struct_1.cond ,'PreLesion')&& day_find > 6)
         taskthreshold = 100000;
     end
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);


mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );

 
 

 %% --------------------------------------------------------------
 case 'JDM25post'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'UniLesion';
descriptor_struct_1.vidtag = 'JDM25_post';

   mocapfilearray=[];
 mocapfiletimes = [];

 for day_find = [7:8]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );


  
 %% --------------------------------------------------------------
 case 'JDM25bipost'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'BiLesion';
descriptor_struct_1.vidtag = 'JDM25_bipost';

   mocapfilearray=[];
 mocapfiletimes = [];

 for day_find = [7:8]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );

   %% --------------------------------------------------------------
 case 'JDM25bipost_task'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'BiLesion';
descriptor_struct_1.vidtag = 'JDM25_bipost_task';

   mocapfilearray=[];
 mocapfiletimes = [];

 for day_find = [7:8,10]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings>10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
  
   case 'JDM25bipost_all'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'BiLesion';
descriptor_struct_1.vidtag = 'JDM25_bipost_all';

   mocapfilearray=[];
 mocapfiletimes = [];

 for day_find = [7:8,10]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);


mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
  
  
  %-------------------------
   case 'JDM25bipost2'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'BiLesion';
descriptor_struct_1.vidtag = 'JDM25_bipost2';

   mocapfilearray=[];
 mocapfiletimes = [];

 for day_find = [10]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );




 %% --------------------------------------------------------------
 case 'JDM25amph'
descriptor_struct_1.tag = 'amph';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM25_amph';

   descriptor_struct_1.day = 8;
   good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds(1:drug_numfiles));
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds(1:drug_numfiles));

 %% --------------------------------------------------------------
 case 'JDM25caff'
descriptor_struct_1.tag = 'caff';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM25_caff';

   descriptor_struct_1.day = 9;
   good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds(1:drug_numfiles));
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds(1:drug_numfiles));



%------------------------
 case 'JDM32pre'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM32_pre';

mocapfilearray = [];
   mocapfiletimes = [];
 good_inds = [];
 
 for day_find = [3:4]
  descriptor_struct_1.day = day_find;
  if ~isempty(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find})
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');
  end
 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );

  %------------------------
 case 'JDM32pre_all'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM32_pre_all';

mocapfilearray = [];
   mocapfiletimes = [];
 good_inds = [];
 
 for day_find = [3:6]
  descriptor_struct_1.day = day_find;
  if ~isempty(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find})
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');
  end
 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );

  
  %------------------------
 case 'JDM32pre2'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM32_pre2';

mocapfilearray = [];
   mocapfiletimes = [];
 good_inds = [];
 
 for day_find = [5:6 ]
  descriptor_struct_1.day = day_find;
  if ~isempty(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find})
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');
  end
 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
  
  
  
  %------------------------
 case 'JDM32pre_task'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM32_pre_task';

mocapfilearray = [];
   mocapfiletimes = [];
 good_inds = [];
 
 for day_find = [3:6 ]
  descriptor_struct_1.day = day_find;
  if ~isempty(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find})
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings>10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');
  end
 end
 
 
 goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );

  %------------------------
 case 'JDM32postuni'
descriptor_struct_1.tag = 'overnight';
descriptor_struct_1.cond = 'UniLesion';
descriptor_struct_1.vidtag = 'JDM32_postuni';

mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [9:10]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
  goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );

  
  
  %------------------------
case 'JDM32postbi'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'BiLesion';
descriptor_struct_1.vidtag = 'JDM32_postbi';

   descriptor_struct_1.day = 8;
   good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds(1:6));
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds(1:6));



%------------------------
case 'JDM32postbi_2'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'BiLesion';
descriptor_struct_1.vidtag = 'JDM32_postbi_2';


mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [8]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
  goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
  

%------------------------
case 'JDM32postbi_3'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'BiLesion';
descriptor_struct_1.vidtag = 'JDM32_postbi_3';


mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [9]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings<10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
  goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
  

  
  %------------------------
case 'JDM32postbi_all'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'BiLesion';
descriptor_struct_1.vidtag = 'JDM32_postbi_all';


mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [9]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
  goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
  
  
  
  %------------------------
case 'JDM32postbi_task'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'BiLesion';
descriptor_struct_1.vidtag = 'JDM32_postbi_task';


mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [8:9]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);
nontaskfiles = find(mocapthreshcrossings>10000);

good_inds = intersect(good_inds,nontaskfiles);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end
 
  goodinds = union(find(datetime(mocapfiletimes).Hour>standard_range_late) , find(datetime(mocapfiletimes).Hour<standard_range_early));
 numsample = min(numel(goodinds),standard_num_files);
 goodinds_sample = randsample(goodinds,numsample);
 
 mocapfilearray = mocapfilearray( goodinds_sample );
  mocapfiletimes = mocapfiletimes( goodinds_sample );
  
  
  
  %------------------------
case 'JDM21'
descriptor_struct_1.tag = 'Rec';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM21_prelesion';


mocapfilearray=[];
 mocapfiletimes = [];
 
 for day_find = [2]
  descriptor_struct_1.day = day_find;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find},descriptor_struct_1.tag)));
mocapfiletimes_temp = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds);
mocapthreshcrossings = mocapfilestruct.(descriptor_struct_1.cond).threshcrossings{day_find}(good_inds);

mocapfilearray = cat(1,mocapfilearray,mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{day_find}(good_inds)');
mocapfiletimes = cat(1,mocapfiletimes,mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{day_find}(good_inds)');

 end

  
  
  %-------------------

 case 'JDM32caff'
descriptor_struct_1.tag = 'caff';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM32_caff';

   descriptor_struct_1.day = 9; %% change this
   good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds(1:drug_numfiles));
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds(1:drug_numfiles));


 case 'JDM32amph'
descriptor_struct_1.tag = 'amph';
descriptor_struct_1.cond = 'PreLesion';
descriptor_struct_1.vidtag = 'JDM32_amph';

   descriptor_struct_1.day = 9; %% change this
   good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds(1:drug_numfiles));
mocapfiletimes = mocapfilestruct.(descriptor_struct_1.cond).mocapdatecreate{descriptor_struct_1.day}(good_inds(1:drug_numfiles));

end

