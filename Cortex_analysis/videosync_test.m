%% load mocap struct
%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';

mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots_videoexamples',filesep);
mkdir(savedirectory);



%% load or create struct
%createmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);
[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_videotest',mocapmasterdirectory);

[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_caff',mocapmasterdirectory);

[mocapstruct_caff] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes,0,0,[],mocapvideodirectory);

[modular_cluster_properties_social] = get_modularclustproperties(mocapstruct_social_2);


[modular_cluster_properties_social] = get_modularclustproperties(mocapstruct_social_2);

%animate_markers_aligned_fullmovie(mocapstruct_social,modular_cluster_properties_social.clustering_inds_agg{2}(1:20:end))
make_cnn_features(mocapvideodirectory)

%% convert video file directories to mp4
% suffix = {'U','L','R'};
% cameratype_ind = 2;
% cameradirectory = strcat(mocapvideodirectory,filesep,'Camera',suffix{cameratype_ind},filesep);
% videofolders = dir(cameradirectory);
% good_directories = find(cat(1, videofolders.isdir)==1);
% 
% for mm = good_directories'
%     if (numel(strfind(videofolders(mm).name,'6')))
%         good_folder = videofolders(mm).name;
%     end
% end
% 
% imageseries_folder = strcat(cameradirectory,good_folder,filesep);
% pre
% 
% %% conver the mkv files to mp4 files
% fprintf('converting mkv files \n')
% camfolder = imageseries_folder;%strcat(imageseries_folder,cameradirectory);
% read_mkv_file(camfolder,camfolder );
% 
% 
% %% load .times file
% %.times file
% times_files = dir(strcat(camfolder,filesep,'*.times'));
% 
% f = fopen(strcat(camfolder,filesep,times_files.name));
% float1 = fread(f,[1,100000000],'uint64');
% frame_number = numel(float1);
% 
% 
% 
% %% synchronize .times file -- associate each frame with a frame
% offset = float1(1);
% video_frames = offset+round((0:size(mocapstruct_social.markers_preproc.HeadF,1)-1)*(1000/300));
% 
% matched_frames =arrayfun(@(x) find(video_frames(x)-float1>0,1,'last'),1:1000000,'UniformOutput', false);
% 
% bad_frames = find(cellfun(@numel,matched_frames) == 0);
% matched_frames_aligned = cat(2,zeros(1,numel(bad_frames)),cell2mat(matched_frames));
% 



  time_subset = find( mocapstruct_social_2.markers_preproc.HeadF(:,3)>160);

time_intersect = intersect(modular_cluster_properties_social.clustering_inds_agg{2},time_subset);


%% choose frames to synchronize,load videos
v = VideoWriter(strcat(savedirectory,'movie example'),'MPEG-4');
open(v)

filename = '\\olveczky.rc.fas.harvard.edu\Jesse\Data\Motionanalysis_captures\Vicon8\20170821\Preprocessed\social_2.mat';
mocapstruct_social = load(filename);

%load in the motion capture struct the field 'matched frames aligned'
%defines the corresponding videoframe for each motion capture frame.
%The video is at 50 Hz MOCAP at 300 Hz, so multiple motion capture
%frames match onto each video frame. In actuality, there is a slight
%mismatch that occurs every 30 minutes (540000 mocap frames) that is
%corrected for below. You can just start with the first 540000 frames,
%where this isn't an issue.
mocapstruct_caff = load('caff.mat');

camerause = 2;
% change the videodirectory
mocapstruct_caff.cameradirectory{camerause} = INSERT VIDEODIRECTORY HERE
%this selects the base frames to start visualization from
base = 1500000;
%this defines a compensation for an odd shift that occurs every 540000
%(30 minutes)
%frames that I am still looking into. The data is synchronized enough
%to use, or you can just start with the first 30 minutes
offset = 10*floor(base./540000);
% this is a visualization of the movies.
M= animate_markers_aligned_fullmovie_syncedvideo(mocapstruct_caff,-offset +(base:10:(base+10000)),...
    mocapstruct_caff.cameradirectory,mocapstruct_caff.matched_frames_aligned{camerause}(base:10:(base+10000)));

writeVideo(v,M)
close(v)
                




%% also save an accelerated video
M_an = animate_markers_aligned_fullmovie(mocapstruct_social,modular_cluster_properties_social.clustering_inds_agg{2}((1:150:end)));

%load data
caffdata = load('E:\Dropbox\mocapdata_for_tim\caff.mat');
%this is just a compatibility issue that will be fixed
caffdata.mocap_struct.rest_frames = 1;
caffdata.mocap_struct.move_frames = 2:size(caffdata.mocap_struct.aligned_mean_position,1);

%get the indicies of the frames when all the markers are tracked
[modular_cluster_properties_caff] = get_modularclustproperties(caffdata.mocap_struct);
%look at cluster 2 - the spine, head and hips
cluster_pick =2;
cluster_pick = 
%add the features to the cluster -- this splices and hi-passes all the data
%together around th
 [modular_cluster_properties_caff2] = get_clustering_features(caffdata.mocap_struct,modular_cluster_properties_caff,cluster_pick) ;

animate_markers_aligned_fullmovie(caffdata.mocap_struct,modular_cluster_properties_caff.clustering_inds_agg{cluster_pick}((1:20:end)));
%

