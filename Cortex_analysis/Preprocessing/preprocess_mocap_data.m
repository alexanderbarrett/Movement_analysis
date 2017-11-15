function [mocap_struct] = preprocess_mocap_data(filepath_array,mocapfilestruct,descriptor_struct,mocapfiletimes,overwrite_flag,overwrite_macro_flag,filepath_array_htr)
%% return a mocapstruct from filepaths
% Jesse Marshall 9/26/2017

[mocap_datecreate,sorted_mocap_serialtimes,filepath_array_sorted ] = sort_mocap_files(filepath_array,mocapfilestruct.mocapdir);

if  nargin ==7
    fprintf('processing HTR files \n')
    [htr_struct] = concat_htr(filepath_array_htr);
end


%% get file directories for data that is saved or not
%% filename for the macrosave
 [save_dir,~,~] = fileparts(filepath_array_sorted{1});
    preproc_save_directory = strrep(save_dir,'Generated_C3D_files','Preprocessed\');
mkdir(preproc_save_directory);
if nargin==7
    macro_save_name = strcat(preproc_save_directory,descriptor_struct.vidtag,'_HTR.mat');
else
macro_save_name = strcat(preproc_save_directory,descriptor_struct.vidtag,'.mat');
end

if (exist(macro_save_name,'file') && ~overwrite_macro_flag)
mocap_struct = load(macro_save_name);
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


%% save the htr info
mocap_struct.limbposition = limbposition;
mocap_struct.rotations = rotations;
mocap_struct.translations_local = translations_local;
mocap_struct.translations_global = translations_global;





%% save this struct
save(save_filename,'mocap_struct','-v7.3')
    end
end
mocap_struct_agg = load_preprocessed_data(filepath_array_sorted);

%% now load from disk and merge


%% level and align the full marker set
%% level and get the aligned
%% make sure the z-axis is pointed up
fprintf('leveling the marker set %f \n')
markers_preproc = level_markers(mocap_struct_agg.markers_preproc);


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
mocap_struct.markercolor = mocapfilestruct.(descriptor_struct.cond).markercolor{descriptor_struct.day};
mocap_struct.links = mocapfilestruct.(descriptor_struct.cond).links{descriptor_struct.day};
mocap_struct.plotdirectory = strcat(mocapfilestruct.mocapdir,(descriptor_struct.cond),'_',num2str(descriptor_struct.day),'_',descriptor_struct.tag);
mkdir(mocap_struct.plotdirectory)


if  nargin ==7
    fprintf('saving HTR files \n')
    htr_fields = fieldnames(htr_struct);
  for  mm = 1:numel(htr_fields)
    mocap_struct.(htr_fields{mm}) = htr_struct.(htr_fields{mm});
  end
 
end

fprintf('saving to: %s \n',macro_save_name);
save(macro_save_name,'mocap_struct','-v7.3')

end
%save(' ',mocap_struct)

end