function [mocap_struct] = preprocess_mocap_data(filepath_array,mocapfilestruct,descriptor_struct)
%% return a mocapstruct from filepaths
% Jesse Marshall 9/26/2017

[mocap_datecreate,sorted_mocap_serialtimes,filepath_array_sorted ] = sort_mocap_files(filepath_array,mocapfilestruct.mocapdir);

% btkCloseAcquisition(acq);
fps = 300;

analog_factor = 1;
analog_fps = fps*analog_factor;

desired_length = 8*10^6;
chunksize = 300*fps;

%filepath_array_sorted_analog = get_analogactive_files(filepath_array_sorted);

[markers,analog,resample_analog,lever_thresholded,filestartpts] = concatenate_andsample_c3d(filepath_array_sorted,fps,analog_fps,...
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


%% get the preprocessed and aligned markrs
[~,markers_preproc_aligned,mean_position,rotation_matrix] = align_hands_elbows(markers_preproc,fps);


fraction_missing = zeros(1,numel(marker_names));
for ll = 1:numel(marker_names)
    missing_times_postpreprocess{ll} = find(markers_preproc.(marker_names{ll})(:,1)==0);
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

% plot(conv(abs_velocity_antialiased(5,:),ones(1,300)./300))

%% get move/not move with a simple threshold
veltrace = (conv(abs_velocity_antialiased(5,:),ones(1,300)./300,'same'));
vel_thresh = 0.01;

frames_move = find(veltrace>vel_thresh);
frames_rest = find(veltrace<=vel_thresh);


mocap_struct.markers_preproc = markers_preproc;
mocap_struct.markers_aligned_preproc = markers_preproc_aligned;
mocap_struct.analog = analog;
mocap_struct.move_frames = frames_move;
mocap_struct.rest_frames = frames_rest;
mocap_struct.bad_frames_agg = bad_frames_agg;
mocap_struct.fraction_missing = fraction_missing;
mocap_struct.resample_analog = resample_analog;
mocap_struct.lever_thresholded = lever_thresholded;
mocap_struct.filestartpts = filestartpts;
mocap_struct.markernames = marker_names;
mocap_struct.fps = fps;
mocap_struct.analog_fps = analog_fps;
mocap_struct.aligned_rotation_matrix = rotation_matrix;
mocap_struct.aligned_mean_position = mean_position;


mocap_struct.markercolor = mocapfilestruct.(descriptor_struct.cond).markercolor{descriptor_struct.day};
mocap_struct.links = mocapfilestruct.(descriptor_struct.cond).links{descriptor_struct.day};
mocap_struct.plotdirectory = strcat(mocapfilestruct.mocapdir,(descriptor_struct.cond),'_',num2str(descriptor_struct.day),'_',descriptor_struct.tag);
mkdir(mocap_struct.plotdirectory)

%save(' ',mocap_struct)

end