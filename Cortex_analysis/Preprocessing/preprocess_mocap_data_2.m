function [mocap_struct] = preprocess_mocap_data(filepath_array,mocapfilestruct,descriptor_struct,mocapfiletimes,overwrite_flag,overwrite_macro_flag,filepath_array_htr,mocapvideodirectory,save_preproc_flag)
%% return a mocapstruct from filepaths
% Jesse Marshall 9/26/2017

[mocap_datecreate,sorted_mocap_serialtimes,filepath_array_sorted ] = sort_mocap_files(filepath_array,mocapfilestruct.mocapdir);

if  nargin ==7
    fprintf('processing HTR files \n')
    [htr_struct] = concat_htr(filepath_array_htr);
end

%% currently only built for a single video directory
% if nargin == 8
%    fprintf('preprocessing the video files \n')
%    
% %% convert video file directories to mp4
% suffix = {'U','L','R'};
% cameratype_ind = 2;
% cameradirectory = strcat(mocapvideodirectory,filesep,'Camera',suffix{cameratype_ind},filesep);
% videofolders = dir(cameradirectory);
% good_directories = find(cat(1, videofolders.isdir)==1);
% camfolder_list = cell(1,numel(good_directories));
% for mm = good_directories'
%     if (numel(strfind(videofolders(mm).name,'6')))
%         good_folder = videofolders(mm).name;
%     end
%     
%     imageseries_folder = strcat(cameradirectory,good_folder,filesep);
% 
% 
% %% conver the mkv files to mp4 files
% fprintf('converting mkv files \n')
% camfolder = imageseries_folder;%strcat(imageseries_folder,cameradirectory);
% read_mkv_file(camfolder,camfolder );
% camfolder_list{mm} = camfolder;
% end
% 
% 
% 
% end




%% get file directories for data that is saved or not
%% filename for the macrosave
 [save_dir,~,~] = fileparts(filepath_array_sorted{1});
    preproc_save_directory = strrep(save_dir,'Generated_C3D_files','Preprocessed\');
mkdir(preproc_save_directory);
if nargin==7
    macro_save_name = strcat(preproc_save_directory,descriptor_struct.Vidtag_savetag,'_HTR.mat');
else
macro_save_name = strcat(preproc_save_directory,descriptor_struct.Vidtag_savetag,'.mat');
end


 %% only save if a new file
    if (exist(macro_save_name,'file'))
[dum,str] = dos(char(strcat('dir',{'   '},macro_save_name)));
dum
c = textscan(str,'%s');
createdate = c{1}{15};
monthval= str2num(createdate(1:2));

if (monthval>1)
save_preproc_flag = 1;
overwrite_macro_flag =1;
fprintf('Overwriting Macro File... \n')
end
end

   fprintf('macros save name %s \n',macro_save_name)
if (exist(macro_save_name,'file') && ~overwrite_macro_flag)
    fprintf('loading from file  \n')
mocap_struct = load(macro_save_name);
if (~isfield(mocap_struct,'markernames'))
mocap_struct = mocap_struct.mocap_struct;
end

[mocap_struct] = assign_modular_annotation_properties(mocap_struct,2);

else
% btkCloseAcquisition(acq);
fps = 300;

analog_factor = 1;
analog_fps = fps*analog_factor;

desired_length = 8*10^6;
chunksize = 300*fps;

%filepath_array_sorted_analog = get_analogactive_files(filepath_array_sorted);

for mm = 1:numel(filepath_array_sorted)
    [save_dir,filename_here,ext_here] = fileparts(filepath_array_sorted{mm});

    preproc_save_directory = strrep(save_dir,'Generated_C3D_files','Preprocessed\');
mkdir(preproc_save_directory);
save_filename= strcat(preproc_save_directory,filename_here,'.mat');

if (exist(save_filename,'file'))
[dum,str] = dos(char(strcat('dir',{'   '},save_filename)));
str
c = textscan(str,'%s');
createdate = c{1}{15};
monthval= str2num(createdate(1:2));
dayval= str2num(createdate(4:5));

if (monthval>1 && dayval > 13)
overwrite_flag = 1;
fprintf('Overwriting ... \n')
end
end

    if (~exist(save_filename,'file') || overwrite_flag)

[markers,analog,resample_analog,lever_thresholded,filestartpts] = concatenate_andsample_c3d(filepath_array_sorted(mm),fps,analog_fps,...
    desired_length,chunksize);


%% get basic information about the dataset
%plot_marker_lever_xcorr(markers,lever_thresholded,analog);


%% touch up the arms so that the elbow to arm is pointed away from the body
%because of elbow and arm swaps, also set to zero if only one of the limbs
%is observed
[markers,markers_aligned,~,~] = align_hands_elbows(markers,fps);



%% Data cleaning
% get frames to analyze based on set conditions
marker_names = fieldnames(markers);
marker_frame_length = size(markers.(marker_names{1}),1);
markers_preproc = markers;
markers_preproc_aligned = markers_aligned;


%% print the 'up times' for markers, and get the missing frames
% get rid of frames where markers are absent
missing_times = cell(1,numel(marker_names));
missing_times_postpreprocess = cell(1,numel(marker_names));

for ll = 1:numel(marker_names)
    missing_times{ll} = find(markers.(marker_names{ll})(:,1)==0);
    fprintf('For marker %s percentage of frames present %e \n',marker_names{ll},...
        100.*numel(missing_times{ll})./marker_frame_length);
end

fprintf('Percentage where two or more are gone %e \n',0)
bad_times = cat(1,missing_times{:});
good_times = setxor(1:marker_frame_length,bad_times);


fakedata = struct();
for fn = fieldnames(markers_preproc)'
    fakedata.(fn{1}) = cell(1,1);
end


%% other interpolation etc.
for ll = 1:numel(marker_names)
    fprintf('starting interpolation for markers, over 30 frames (100 ms) maximum %f \n',ll)
    % for jkj = 1:3 %need to do over all x y z simul and add spike correction
    [temp,fake_frames] = markershortinterp((markers_preproc.(marker_names{ll})),30,5);
    fakedata.(marker_names{ll}){1} = fake_frames;
    % if numel(temp)
    markers_preproc.(marker_names{ll}) = temp;
    %  end
    clear fake_frames
end




% for ll = 1:numel(marker_names)
%     fprintf('interpolating aligned markers \n')
%     [temp,~] =markershortinterp((markers_preproc.(marker_names{ll})),100,5);
%     % if numel(temp)
%     markers_preproc_aligned.(marker_names{ll}) = temp;
%     %  end
% end

%% median filter the data to remove spikes
fprintf('median filtering %f \n')

for ll = 1:numel(marker_names)
    markers_preproc.(marker_names{ll}) = medfilt2(markers_preproc.(marker_names{ll}),[3,1]);
end


%print('-dpng',strcat(plotdirectory,'MarkerTracking.png'))
%set(get(gca,'XTickLabel'),'Rotation',-45);



%% Transform the data into relevant features for clustering and visualization

%'big-data' features
% get relative marker positions to one another (x,y,z)
num_markers = numel(markers);
marker_velocity = zeros(num_markers,marker_frame_length,4);
marker_position = zeros(num_markers,marker_frame_length,3);
abs_velocity_antialiased = zeros(num_markers,marker_frame_length);


dH = designfilt('lowpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 60/(fps/2), ...
    'DesignMethod', 'butter');
[f1,f2] = tf(dH);

%delta_markers_reshaped = [];
fprintf('getting velocities \n')
for ll = 1:numel(marker_names)
    marker_position(ll,:,1:3) = markers_preproc.(marker_names{ll});
    for mk = 1:3
    end   
    marker_velocity(ll,2:(end),1:3) = diff(markers_preproc.(marker_names{ll}));
    marker_velocity(ll,1,1:3) = marker_velocity(ll,2,1:3);
    marker_velocity(ll,:,4) = sqrt(sum((squeeze( marker_velocity(ll,:,1:3))).^2,2));
    abs_velocity_antialiased(ll,:) =  filtfilt(f1,f2, marker_velocity(ll,:,4));

    for lk = (ll+1):num_markers
        distance_here =   (markers_preproc.(marker_names{ll})-markers_preproc.(marker_names{lk}));
    end
end


%get aggregate feature matrix
%% simple bad frame detector
fprintf('finding bad frames \n')
bad_frames_agg = getbadframes(marker_velocity,marker_position,fps);
clear marker_position



%% get rest/move
%% get move/not move with a simple threshold
veltrace = (conv(abs_velocity_antialiased(5,:),ones(1,300)./300,'same'));
vel_thresh = 0.008;
vel_thresh_fast = 0.05;

frames_move = find(veltrace>vel_thresh);
frames_rest = find(veltrace<=vel_thresh);


frames_move_fast = find(veltrace>vel_thresh_fast);
frames_rest_fast = find(veltrace<=vel_thresh_fast);


move_frames = zeros(1,numel(veltrace));
move_frames(veltrace>vel_thresh) = 1;
move_frames = conv(move_frames,ones(1,300)./300,'same');

frames_near_move = find(move_frames >0);
frames_near_rest = setxor(1:numel(veltrace),frames_move);


%% move/not move -- can change
mocap_struct.move_frames = frames_move;
mocap_struct.rest_frames = frames_rest;
mocap_struct.move_frames_fast = frames_move_fast;
mocap_struct.rest_frames_fast = frames_rest_fast;
mocap_struct.move_near_frames = frames_near_move;
mocap_struct.rest_near_frames = frames_near_rest;
% plot(conv(abs_velocity_antialiased(5,:),ones(1,300)./300))

%% fill the struct fields
mocap_struct.analog = analog;
mocap_struct.bad_frames_agg = bad_frames_agg;
mocap_struct.resample_analog = resample_analog;
mocap_struct.lever_thresholded = lever_thresholded;
mocap_struct.markers_preproc = markers_preproc;
mocap_struct.filestartpts = filestartpts;

mocap_struct.filenames = filepath_array_sorted(mm);

%% save the htr info

if  nargin ==7
    fprintf('saving HTR files \n')
    htr_fields = fieldnames(htr_struct);
  for  mm = 1:numel(htr_fields)
    mocap_struct.(htr_fields{mm}) = htr_struct.(htr_fields{mm});
  end
 
end

%% do the video file analysis

[dir_base_file,fname_here] = fileparts(filepath_array_sorted{mm});
%find which condition the file is from
day_here = nan;
condhere = [];

days_to_search = descriptor_struct.Days;
if (ischar(descriptor_struct.Days))
days_to_search = str2num( days_to_search);
end
for daysearch =days_to_search
condcandidates = mocapfilestruct.(descriptor_struct.Condition).day_conds{daysearch};
for mm = condcandidates
    
    mocapfileshere = cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct.Condition).mocapfiles{daysearch},mm));
mocapfile_exact = cellfun(@numel,(strfind(mocapfilestruct.(descriptor_struct.Condition).mocapfiles{daysearch},fname_here)));
[~,mocapfile_number,~] = intersect(find(mocapfileshere),find(mocapfile_exact));
    
    
    if numel(mocapfile_number)
        day_here = daysearch;
        condhere = mm;
    end
end
end

if ~isnan(day_here)
  %% do the video part here
%videofiles
dir_base = ((strrep(dir_base_file,'Generated_C3D_files','')));
file_folders = dir(dir_base);
directories = find(cat(1,file_folders.isdir)==1);
good_folder = [];
for mm = directories'
if (numel(strfind(file_folders(mm).name,condhere)))
    good_folder = file_folders(mm).name;
end
end
mocapvideodirectory =  strcat(dir_base,good_folder);
  %------------------------

   fprintf('preprocessing the video files \n')
   
%% convert video file directories to mp4
suffix = {'U','L','R'};
camfolder_list = cell(1,numel(suffix));

for cameratype_ind = 1:numel(suffix)
cameradirectory = strcat(mocapvideodirectory,filesep,'Camera',suffix{cameratype_ind},filesep);
if exist(cameradirectory,'dir')
videofolders = dir(cameradirectory);
for mm =1:numel(videofolders)

%good_directories = find(cat(1, videofolders.isdir)==1);
    if (numel(strfind(videofolders(mm).name,'6')))
        good_folder = videofolders(mm).name;
    
    imageseries_folder = strcat(cameradirectory,good_folder,filesep);


%% conver the mkv files to mp4 files
fprintf('converting mkv files \n')
camfolder = imageseries_folder;%strcat(imageseries_folder,cameradirectory);
%convert mkv to mp4 (slow)
%read_mkv_file(camfolder,camfolder );
camfolder_list{cameratype_ind} = camfolder;
    end

end
end
end

%% align the frames
%first get the mocapfiles for this condition
mocapfileshere = cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct.Condition).mocapfiles{day_here},condhere));
mocapfile_exact = cellfun(@numel,(strfind(mocapfilestruct.(descriptor_struct.Condition).mocapfiles{day_here},fname_here)));
[~,mocapfile_number,~] = intersect(find(mocapfileshere),find(mocapfile_exact));
%then the cam directory
mocap_struct.cameradirectory = cell(1,numel(camfolder_list));
if numel(camfolder_list) == 0
    fprintf('no camera folder found \n')
       mocap_struct.matched_frames_aligned{mm} = {[],[],[]}; 
end

    for mm=1:numel(camfolder_list)
        camfolder = camfolder_list{mm};
        if numel(camfolder)
   fprintf('saving video properties \n')
mocap_struct.cameradirectory{mm} =  camfolder;

%% load .times file
times_files = dir(strcat(camfolder,filesep,'*.times'));

f = fopen(strcat(camfolder,filesep,times_files.name));
fprintf('loaded times file %f \n',f);

if (f)
float1 = fread(f,[1,100000000],'uint64');
frame_number = numel(float1);

%% synchronize .times file -- associate each frame with a frame
if numel(float1)
offset_video = float1(1);
float1 = float1-offset_video;

%get offset of the filenumber
%mocapfile_number, the first is standard number of frames, the second is
%ann odd shift that occurs with each file number
offset_mocap = (mocapfile_number-1)*540000*(1000/300)+780*(mocapfile_number-1)*(1000/300);
video_frames = bsxfun(@plus,round((0:size(mocap_struct.markers_preproc.HeadF,1)-1)*(1000/300)),offset_mocap);

matched_frames_last =arrayfun(@(x) find(video_frames(x)-float1>0,1,'last'),numel(video_frames),'UniformOutput', false);
matched_frames_first =arrayfun(@(x) find(video_frames(x)-float1>0,1,'last'),1,'UniformOutput', false);
if numel(matched_frames_first{1}) ==0
    matched_frames_first{1}  = 1;
end
reduced_float1 = float1(matched_frames_first{1}:matched_frames_last{1});
matched_frames =arrayfun(@(x) (matched_frames_first{1}-1)+find(video_frames(x)-reduced_float1>0,1,'last'),1:numel(video_frames),'UniformOutput', false);



bad_frames = find(cellfun(@numel,matched_frames) == 0);
matched_frames_aligned = cat(2,zeros(1,numel(bad_frames)),cell2mat(matched_frames));

% get the matched frames here
mocap_struct.matched_frames_aligned{mm} = matched_frames_aligned;
else
    mocap_struct.matched_frames_aligned{mm} = [];
end
        else
            fprintf('no camera folders found \n')
             mocap_struct.matched_frames_aligned{mm} = [];
end
else
     mocap_struct.matched_frames_aligned{mm} = [];
end
%% load in all the video files and save to the structure
%mocap_struct.video_files = 
    end
else
    fprintf('no camera folders found \n')
             mocap_struct.matched_frames_aligned = {[],[],[]};
                          mocap_struct.cameradirectory = {[],[],[]};

end

        fprintf(' adding outputs \n')
       [~,markers_preproc_aligned,mean_position,rotation_matrix] = align_hands_elbows(mocap_struct.markers_preproc,fps);
mocap_struct.markercolor = mocapfilestruct.(descriptor_struct.Condition).markercolor{1}; 
     mocap_struct.links = mocapfilestruct.(descriptor_struct.Condition).links{1};
mocap_struct.markernames = fieldnames(mocap_struct.markers_preproc ); 
  mocap_struct.markers_aligned_preproc = markers_preproc_aligned;
  mocap_struct.fps = fps;
mocap_struct.analog_fps = analog_fps;
    mocap_struct.mocapfiletimes = mocapfiletimes;


%% save this struct
save(save_filename,'mocap_struct','-v7.3')
    else
       load(save_filename)
       if ~isfield(mocap_struct,'markers_aligned_preproc')
           fprintf(' adding outputs \n')
       [~,markers_preproc_aligned,mean_position,rotation_matrix] = align_hands_elbows(mocap_struct.markers_preproc,fps);
mocap_struct.markercolor = mocapfilestruct.(descriptor_struct.Condition).markercolor{1}; 
     mocap_struct.links = mocapfilestruct.(descriptor_struct.Condition).links{1};
mocap_struct.markernames = fieldnames(mocap_struct.markers_preproc ); 
  mocap_struct.markers_aligned_preproc = markers_preproc_aligned;
  mocap_struct.fps = fps;
mocap_struct.analog_fps = analog_fps;
mocap_struct.mocapfiletimes = mocapfiletimes;

       save(save_filename,'mocap_struct')
       end
    end
end
mocap_struct_agg = load_preprocessed_data(filepath_array_sorted);

%% add fields left out by a previous configuration

for ll = 1:size(mocap_struct_agg.cameradirectory,2)
    fprintf('converting mkv files')
    read_mkv_file(mocap_struct.cameradirectory{1,ll},mocap_struct.cameradirectory{1,ll});
end

%% now load from disk and merge


%% level and align the full marker set
%% level and get the aligned
%% make sure the z-axis is pointed up
fprintf('leveling the marker set %f \n')
if numel(fieldnames( mocap_struct_agg.markers_preproc)) == 20
markers_preproc = level_markers(mocap_struct_agg.markers_preproc);
else
markers_preproc = mocap_struct_agg.markers_preproc;
end
%% get the preprocessed and aligned markrs
[~,markers_preproc_aligned,mean_position,rotation_matrix] = align_hands_elbows(mocap_struct_agg.markers_preproc,fps);
marker_names = fieldnames(mocap_struct_agg.markers_preproc );
marker_frame_length = size(mocap_struct_agg.markers_preproc.(marker_names{1}),1);


%% report the fraction missing
fraction_missing = zeros(1,numel(marker_names));
for ll = 1:numel(marker_names)
    missing_times_postpreprocess{ll} = find(mocap_struct_agg.markers_preproc.(marker_names{ll})(:,1)==0);
    fprintf('For marker %s percentage of missing frames %e \n',...
        marker_names{ll},numel(missing_times_postpreprocess{ll})./marker_frame_length);
    fraction_missing(ll) = numel(missing_times_postpreprocess{ll})./marker_frame_length;
end

figure(555)
bar(cellfun(@numel,missing_times_postpreprocess)./marker_frame_length)
box off
set(gca,'XTick',1:numel(marker_names),'XTickLabels',marker_names)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Fraction of time tracked')
ylim([0 0.5])
xlim([0 numel(marker_names)+1])




%% save the matched video frames
% Now done for individual files
% 
% if  nargin ==8
%     for mm=1:1%numel(camfolder_list)
%         camfolder = camfolder_list{mm};
%    fprintf('saving video properties \n')
% mocap_struct.cameradirectory = cat(2,mocap_struct.cameradirectory,camfolder);
% 
% %% load .times file
% times_files = dir(strcat(camfolder,filesep,'*.times'));
% 
% f = fopen(strcat(camfolder,filesep,times_files.name));
% float1 = fread(f,[1,100000000],'uint64');
% frame_number = numel(float1);
% 
% %% synchronize .times file -- associate each frame with a frame
% offset = float1(1);
% video_frames = offset+round((0:size(markers_preproc_aligned.HeadF,1)-1)*(1000/300));
% 
% matched_frames =arrayfun(@(x) find(video_frames(x)-float1>0,1,'last'),1:size(mean_position,1),'UniformOutput', false);
% 
% bad_frames = find(cellfun(@numel,matched_frames) == 0);
% matched_frames_aligned = cat(2,zeros(1,numel(bad_frames)),cell2mat(matched_frames));
% 
% mocap_struct.matched_frames_aligned = matched_frames_aligned;
% 
% %% load in all the video files and save to the structure
% %mocap_struct.video_files = 
%     end
% end





%% get move/not move with a simple threshold
% veltrace = (conv(abs_velocity_antialiased(5,:),ones(1,300)./300,'same'));
% vel_thresh = 0.008;
% vel_thresh_fast = 0.05;
% 
% frames_move = find(veltrace>vel_thresh);
% frames_rest = find(veltrace<=vel_thresh);
% 
% 
% frames_move_fast = find(veltrace>vel_thresh_fast);
% frames_rest_fast = find(veltrace<=vel_thresh_fast);
% 
% 
% move_frames = zeros(1,numel(veltrace));
% move_frames(veltrace>vel_thresh) = 1;
% move_frames = conv(move_frames,ones(1,300)./300,'same');
% 
% frames_near_move = find(move_frames >0);
% frames_near_rest = setxor(1:numel(veltrace),frames_move);
% 
% 
% %% move/not move -- can change
mocap_struct.move_frames = mocap_struct_agg.move_frames;
mocap_struct.rest_frames = mocap_struct_agg.rest_frames;
mocap_struct.move_frames_fast = mocap_struct_agg.move_frames_fast;
mocap_struct.rest_frames_fast = mocap_struct_agg.rest_frames_fast;
mocap_struct.move_near_frames = mocap_struct_agg.move_near_frames;
mocap_struct.rest_near_frames = mocap_struct_agg.rest_near_frames;
mocap_struct.fraction_missing = fraction_missing;

%% doing after for leveling
mocap_struct.markers_aligned_preproc = markers_preproc_aligned;
mocap_struct.markers_preproc = markers_preproc;

mocap_struct.aligned_rotation_matrix = rotation_matrix;
mocap_struct.aligned_mean_position = mean_position;

%% merged prop
mocap_struct.analog = mocap_struct_agg.analog;
mocap_struct.bad_frames_agg = mocap_struct_agg.bad_frames_agg;
mocap_struct.resample_analog = mocap_struct_agg.resample_analog;
mocap_struct.lever_thresholded = mocap_struct_agg.lever_thresholded;
mocap_struct.filestartpts = mocap_struct_agg.filestartpts;


%% global properties -- same for all files

mocap_struct.markernames = marker_names;
mocap_struct.fps = fps;
mocap_struct.analog_fps = analog_fps;
mocap_struct.mocapfiletimes = mocapfiletimes;
mocap_struct.markercolor = mocapfilestruct.(descriptor_struct.Condition).markercolor{1};
mocap_struct.links = mocapfilestruct.(descriptor_struct.Condition).links{1};
mocap_struct.plotdirectory = strcat(mocapfilestruct.mocapdir,(descriptor_struct.Condition),'_',num2str(1),'_',descriptor_struct.Nametag);

if numel(fieldnames( mocap_struct_agg.markers_preproc)) == 20
[mocap_struct] = assign_modular_annotation_properties(mocap_struct,2);
end

mkdir(mocap_struct.plotdirectory)

if  nargin ==7
    fprintf('saving HTR files \n')
    htr_fields = fieldnames(htr_struct);
  for  mm = 1:numel(htr_fields)
    mocap_struct.(htr_fields{mm}) = htr_struct.(htr_fields{mm});
  end
 
end
if (nargin < 9 || save_preproc_flag)
   
    
    
    
fprintf('saving to: %s \n',macro_save_name);
save(macro_save_name,'mocap_struct','-v7.3')
end
end
%save(' ',mocap_struct)

end