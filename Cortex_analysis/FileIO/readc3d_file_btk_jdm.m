% read_c3d_file_better

%% Dependencies

% download btk from: https://code.google.com/archive/p/b-tk/downloads
% unzip, and addpath to matlab

% tutorial
% http://biomechanical-toolkit.github.io/docs/Wrapping/Matlab/_tutorial.html

% function list
% http://biomechanical-toolkit.github.io/docs/Wrapping/Matlab/annotated.html

%% read c3df

%filepath = 'E:\Bence\Data\MOCAP\Rat_video_datasets\20170220\Mocap\TOYDOG_videosync_long1.c3d';
fileselect = 7;
%filepath = 'E:\Bence\Data\Motionanalysis_captures\20170227\Vicon6_task1\Vicon6_task1.c3d';
%filepath = 'E:\Bence\Data\Motionanalysis_captures\20170227\Vicon6_novelobj1\Vicon6_novelobj1.c3d';
%filepath = 'E:\Bence\Data\Motionanalysis_captures\20170228\Vicon6_longrun2\Vicon6_longrun2_postprocessed.c3d';

%filepath = 'E:\Bence\Data\Motionanalysis_captures\20170301\Vicon6_postanesthesia_ketamine1\Vicon6_postanesthesia_ketamine1_postprocessing_firstpass.c3d' ;

%filepath = 'E:\Bence\Data\Motionanalysis_captures\20170222\Vicon12_example_package3-GS-1\Vicon12_example_package3-GS-1_JDM.c3d';
%
%filepath = 'E:\Bence\Data\Motionanalysis_captures\20170228\Vicon6_longrun2\Vicon6_longrun2_postprocessed_badgrooming_joined.c3d';
%
%
%




%% Vicon6 with elbows
filepath = 'E:\Bence\Data\Motionanalysis_captures\20170404\Vicon6_task_long2\Vicon6_task_long2_postp.c3d';

%videofeatures_files = {'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU\videofeaturesU.mat',...
%   'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL\videofeaturesL.mat'};

videofeatures_directory = {'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU',...
    'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL'};
feature_tag =  'videofeatures_convlayer_pca*';
savedirectory = 'E:\Bence\Data\MOCAP\regression_animations\';
save_tag = 'Vicon6_elbows_task_100clusters';
    clustering_ind = [65000:410000];%setxor(1:frame_length_here,intersect(1:frame_length_here,exclude_frames_clustering));

    
%% Vicon6 with one elbow
filepath = 'E:\Bence\Data\Motionanalysis_captures\20170406\Vicon6_rat_kinect_long1\Vicon6_rat_kinect_long1_firstpostpr_interp.c3d';

%videofeatures_files = {'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU\videofeaturesU.mat',...
%   'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL\videofeaturesL.mat'};

videofeatures_directory = {'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU',...
    'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL'};
feature_tag =  'videofeatures_convlayer_pca*';
savedirectory = 'E:\Bence\Data\MOCAP\regression_animations\';
save_tag = 'Vicon6_kinectlong1';
    clustering_ind = [1:440000];%setxor(1:frame_length_here,intersect(1:frame_length_here,exclude_frames_clustering));


        
%% Vicon6 with hands
%filepath = 'E:\Bence\Data\Motionanalysis_captures\20170411\Vicon6_openfield_task_long1\Vicon6_openfield_task_long1_postp_badhand_forexport.c3d';
%filepath = 'E:\Bence\Data\Motionanalysis_captures\20170411\Vicon6_openfield_task_long1\Vicon6_openfield_task_long1_postp_badhand_forexport.c3d';
filepath = 'E:\Bence\Data\Motionanalysis_captures\20170411\Vicon6_longruntask21\Vicon6_longruntask21_shoddyinterp.c3d';

filepath_array = { 'E:\Bence\Data\Motionanalysis_captures\20170411\Vicon6_longruntask21\Vicon6_longruntask21_shoddyinterp.c3d',...
    'E:\Bence\Data\Motionanalysis_captures\20170411\Vicon6_openfield_task_long1\Vicon6_openfield_task_long1_postp_badhand_forexport.c3d',...
    'E:\Bence\Data\Motionanalysis_captures\20170411\Vicon6_longruntask_notask1\Vicon6_longruntask_notask1_postp.c3d'};

%videofeatures_files = {'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU\videofeaturesU.mat',...
%   'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL\videofeaturesL.mat'};

videofeatures_directory = {'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU',...
    'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL'};
feature_tag =  'videofeatures_convlayer_pca*';
savedirectory = 'E:\Bence\Data\MOCAP\regression_animations\';
save_tag = 'Vicon6_hands_aggregated';
    clustering_ind = [1:2400000];%setxor(1:frame_length_here,intersect(1:frame_length_here,exclude_frames_clustering));

    
    
    
    
    %% Vicon3 with one hand
    
mocapdir = 'E:\Bence\Data\Motionanalysis_captures\Vicon3\';
mocapnames = {'20170710\Generated_C3D_files\'};

filenames_here = dir(strcat(mocapdir,mocapnames{1},'Vicon3_recording_636353114036170264*.c3d'));
filepath_array = [];
filepath_ind = 1;
for kk = 1:numel(filenames_here)
    if (~numel(strfind(filenames_here(kk).name,'nolj')))
        filepath_array{filepath_ind} = strcat(mocapdir,mocapnames{1},filenames_here(kk).name);
        filepath_ind = filepath_ind+1;
    end
end
savedirectory = 'E:\Bence\Data\MOCAP\regression_animations\';
save_tag = 'Vicon3_subset_ptii';
    %clustering_ind = [1:12000000];%setxor(1:frame_length_here,intersect(1:frame_length_here,exclude_frames_clustering));

    
    %% Vicon3 with 18 markers
    
mocapdir = 'E:\Bence\Data\Motionanalysis_captures\Vicon3\';
mocapnames = {'20170721\Generated_C3D_files\'};

filenames_here = dir(strcat(mocapdir,mocapnames{1},'Vicon3_recording_amphetamine*.c3d'));
filepath_array = [];
filepath_ind = 1;
for kk = 1:numel(filenames_here)
    if (~numel(strfind(filenames_here(kk).name,'nolj')))
        filepath_array{filepath_ind} = strcat(mocapdir,mocapnames{1},filenames_here(kk).name);
        filepath_ind = filepath_ind+1;
    end
end
savedirectory = 'E:\Bence\Data\MOCAP\regression_animations\';
save_tag = 'Vicon3_18marker_amph3';
    %clustering_ind = [1:12000000];%setxor(1:frame_length_here,intersect(1:frame_length_here,exclude_frames_clustering));
    
    
    
%     
%     
% filepath_array = { 'E:\Bence\Data\Motionanalysis_captures\20170420\Vicon6_longrun1-Vicon6Elbow.c3d',...
%     'E:\Bence\Data\Motionanalysis_captures\20170420\Vicon6_longrun2-Vicon6Elbow.c3d',...
%     'E:\Bence\Data\Motionanalysis_captures\20170420\Vicon6_longrun4-Vicon6Elbow.c3d',...
%     'E:\Bence\Data\Motionanalysis_captures\20170420\Vicon6_longrun5-Vicon6Elbow.c3d',...
%     'E:\Bence\Data\Motionanalysis_captures\20170420\Vicon6_longrun6-Vicon6Elbow.c3d'};
% % 




%     % %% Vicon 6 no hands, no elbows
% filepath = 'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\Vicon6_freemoveexplore1.c3d';
% %videofeatures_files = {'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU\videofeaturesU.mat',...
% %   'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL\videofeaturesL.mat'};
% videofeatures_directory = {'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU',...
%     'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL'};
% feature_tag =  'videofeatures_convlayer_pca*';
% savedirectory = 'E:\Bence\Data\MOCAP\regression_animations\';
% save_tag = 'Vicon6_openfield_noelbows';
%     clustering_ind = [1:400000];%setxor(1:frame_length_here,intersect(1:frame_length_here,exclude_frames_clustering));
%     
%     
%     
%     %% Vicon 6 ketamine
% filepath = 'E:\Bence\Data\Motionanalysis_captures\20170301\Vicon6_postanesthesia_ketamine1\Vicon6_postanesthesia_ketamine1_postprocessing_firstpass.c3d';
% %videofeatures_files = {'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU\videofeaturesU.mat',...
% %   'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL\videofeaturesL.mat'};
% videofeatures_directory = {'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraU',...
%     'E:\Bence\Data\Motionanalysis_captures\20170224\Vicon6_freemoveexplore1\CameraL'};
% feature_tag =  'videofeatures_convlayer_pca*';
% savedirectory = 'E:\Bence\Data\MOCAP\regression_animations\';
% save_tag = 'Vicon6_openfield_kx';
%     clustering_ind = [1:1200000];%setxor(1:frame_length_here,intersect(1:frame_length_here,exclude_frames_clustering));
% 
%     
    
   % filepath = 'E:\Bence\Data\Motionanalysis_captures\20170301\Vicon6_postanesthesia_ketamine1\Vicon6_postanesthesia_ketamine1_postprocessing.c3d';

    
%% TOYDOG RD 1
% filepath = 'E:\Bence\Data\Motionanalysis_captures\20170224\Toydog_long1\Toydog_long1.c3d';
% %videofeatures_files = {'E:\Bence\Data\Motionanalysis_captures\20170224\Toydog_long1\CameraL\videofeaturesL.mat',...
% %    'E:\Bence\Data\Motionanalysis_captures\20170224\Toydog_long1\CameraU\videofeaturesU.mat'};
%
% videofeatures_directory = {'E:\Bence\Data\Motionanalysis_captures\20170224\Toydog_long1\CameraL',...
%     'E:\Bence\Data\Motionanalysis_captures\20170224\Toydog_long1\CameraU'};
% %feature_tag = 'videofeatures_outputlayer*';
% feature_tag = 'videofeatures_convlayer_pca*';
%save_tag = 'Toy_dog_1';

% TOY DOG RD 2
% filepath = 'E:\Bence\Data\Motionanalysis_captures\20170309\toydog_good2\toydog_good2_postprocessed.c3d';
%
% videofeatures_directory = {'E:\Bence\Data\Motionanalysis_captures\20170309\toydog_good2\CameraL',...
%     'E:\Bence\Data\Motionanalysis_captures\20170309\toydog_good2\CameraU','E:\Bence\Data\Motionanalysis_captures\20170309\toydog_good2\CameraR'};
%
% feature_tag =  'videofeatures_convlayer_pca*';
% save_tag = 'Toy_dog_2_decisiontree';
% savedirectory = 'E:\Bence\Data\MOCAP\regression_animations\';


% feature_tag = 'videofeatures_outputlayer*';




toy_dog = 0;
toydog = 0;
do_tsne = 0;
rat_hands = 0;
rat_elbows = 0;
rat_hands_offset = 0;
rat_nohands_noelbows = 0;
vicon3_singlehand = 0;
vicon3_18marker = 1;


num_pcs = 40;
do_video_analysis = 0;
num_clusters =200;
fps = 300;
clustering_window = floor(fps);
clustering_overlap = floor(fps./2);
do_conv_features = 0;
time_fraction_clustering = 1.0;
exclude_frames_clustering = 0;%cat(2,175000:260000,135000:150000,618000:624000);
%params for toydog: 20 clusters
markercolor = {'b','b','r','r','r','m','g','g','k','k','k','k','k'};
links = {[1,2],[2,3],[3,4],[5,6],[4,5],[5,7],[5,8],[4,6],[3,6],[7,8]};



if (vicon3_18marker)
    %elbow is black
    markercolor = {'b','b','b',...
        'r','r','r',...
        'm','m',...
        'g','g',... %hips
        'y','y','y',... %L arm
        'w','w','w',... %R arm
        'g','g',...
        'b','b','k','k','k',};
    links = {[1,2],[2,3],[1,3],...
        [2,6],[1,6],[3,6],... %head to spine
        [6,4],[5,4],...
        [6,7],[7,8],[4,8],[4,7],[5,7],...
        [5,9],[5,10],...
        [11,12],[6,13],[6,14],[11,13],[12,13],...
        [14,15],[14,16],[15,16],...
        [9,18],[10,17]};
end



if (vicon3_singlehand)
    markercolor = {'b','b','b','r','r','r','m','m','g','g','y','y','k','k','k'};
    links = {[1,2],[2,3],[1,3],...
        [3,6],[6,5],[5,4],...
        [6,7],[7,8],[6,8],[4,8],[4,7],[5,7],[5,8],...
        [5,9],[5,10],...
        [6,11],[11,12]};
end

if (rat_hands)
    markercolor = {'b','b','r','r','r','m','g','g','y','y','y','y'};
    links = {[1,2],[2,3],[3,4],[4,6],[5,7],[5,6],[4,5],[5,7],[5,8],[4,6],[3,6],[7,8],[3,9],[3,11],[9,10],[11,12]};
end

if (rat_hands_offset)
    markercolor = {'b','b','r','m','r','m','r','g','g','y','y','y','y'};
    links = {[1,2],[2,3],[3,4],[4,6],[5,7],[8,9],[5,6],[3,5],[4,5],[6,7],[7,8],[7,9],[3,11],[3,10],[10,12],[11,13]};
end


if (rat_elbows)
    markercolor = {'b','b','r','m','r','m','r','g','g','y','y','k','k'};
    links = {[1,2],[2,3],[3,5],[5,7],[3,4],[4,5],[5,6],[6,7],[7,8],[7,9],[8,9],...
        [3,10],[3,11]};
end
        

if (rat_nohands_noelbows)
 markercolor = {'b','b','r','r','r','m','g','g','k','k','k','k','k'};
links = {[1,2],[2,3],[3,4],[5,6],[4,5],[5,7],[5,8],[4,6],[3,6],[7,8]};
        end

if (toy_dog)
    markercolor = {'b','b','r','m','r','g','g','y','y'};
    links = {[1,2],[2,3],[3,5],[3,4],[5,4],[3,6],[6,7],[3,9],[8,9]};
end




%filepath = 'E:\Bence\Data\Motionanalysis_captures\20170301\Vicon6_postanesthesia_ketamine1\Vicon6_postanesthesia_ketamine1_postprocessing.c3d';

%filepath = '\\OLVECZKYLABSERV\Files\Kevin\20161208\Vicon2_mid_3_2.c3d';
% acq = btkReadAcquisition(filepath);
%
% markers = btkGetMarkers(acq);
% analog = btkGetAnalogs(acq);
%
% btkCloseAcquisition(acq);
analog_factor = 1;
analog_fps = fps*analog_factor;

%[markers,analog,resample_analog,lever_thresholded] = readc3d_jdm(filepath,fps,analog_fps);
[markers,analog,resample_analog,lever_thresholded] = concatenate_c3d(filepath_array,fps,analog_fps);


% %% get frame rate
% fid = fopen(filepath, 'r');
% fseek(fid,0,'cof');
% frame_rate=fread(fid,1,'float32');
% fclose(fid);
% fps = 245;
%
%
% %% Process the analog data
% analog_fps = fps*20;
% w0 = 60/(analog_fps/2); %filter fundamental frequency
% bw = w0/50; %Q=35
% [b,a] = iirnotch(w0,bw);
% filtered_lever = filtfilt(b,a,analog.LEVER);
%
% lever_threshold = 0.6;
%
% lever_thresholded = zeros(1,numel(filtered_lever));
% lever_thresholded(filtered_lever>lever_threshold) = 1;
% lever_thresholded(filtered_lever<-lever_threshold) = -1;
%
% figure(111)
% plot(filtered_lever)
% hold on
% plot(0.6*ones(1,numel(filtered_lever)),'r')
% plot(analog.LICK_SENSOR,'g')
% plot(lever_thresholded,'k')
% plot(analog.SPEAKER,'y')
% hold off
%
%
% resample_analog = resample([lever_thresholded' analog.SPEAKER analog.WATER_VALVE analog.LICK_SENSOR analog.GROUND],1,20);
% resample_analog = bsxfun(@minus,resample_analog,mean(resample_analog,1));
%filter the 60 hz noise in the lever


%% get basic information about the dataset
marker_names = fieldnames(markers);
num_markers = numel(marker_names);



for ll = 1:numel(marker_names)
    %markers.(marker_names{ll}) = markers.(marker_names{ll})(1:12000000,:);
end
clustering_ind = 1:size(markers.(marker_names{1})  ,1);
fprintf('number of frames %f \n',size(markers.(marker_names{1}),1));

%% Data cleaning
% get frames to analyze based on set conditions
marker_frame_length = size(markers.(marker_names{1}),1);
markers_preproc = markers;


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
    [temp,fake_frames] =markershortinterp((markers_preproc.(marker_names{ll})),30,5);
    fakedata.(marker_names{ll}){1} = fake_frames;
   % if numel(temp)
    markers_preproc.(marker_names{ll}) = temp;
  %  end
    clear fake_frames
  %  end
end

%% median filter the data to remove spikes
    fprintf('median filtering %f \n')

for ll = 1:numel(marker_names)
    markers_preproc.(marker_names{ll}) = medfilt2(markers_preproc.(marker_names{ll}),[3,1]);
end


for ll = 1:numel(marker_names)
    missing_times_postpreprocess{ll} = find(markers_preproc.(marker_names{ll})(:,1)==0);
    fprintf('For marker %s percentage of missing frames %e \n',...
        marker_names{ll},numel(missing_times_postpreprocess{ll})./marker_frame_length);
end



%% Transform the data into relevant features for clustering and visualization
% biologically inspired features


%'big-data' features
% get relative marker positions to one another (x,y,z)
delta_markers = zeros(num_markers,num_markers,marker_frame_length,4);
marker_velocity = zeros(num_markers,marker_frame_length,4);
marker_position = zeros(num_markers,marker_frame_length,3);

%delta_markers_reshaped = [];
fprintf('getting velocities \n')
for ll = 1:numel(marker_names)
    marker_position(ll,:,1:3) = markers_preproc.(marker_names{ll});
    marker_velocity(ll,2:(end),1:3) = diff(markers_preproc.(marker_names{ll}));
    marker_velocity(ll,1,1:3) = marker_velocity(ll,2,1:3);
    marker_velocity(ll,:,4) = sqrt(sum((squeeze( marker_velocity(ll,:,1:3))).^2,2));
    
    for lk = (ll+1):num_markers
        delta_markers(ll,lk,:,1:3) =   markers_preproc.(marker_names{ll})-markers_preproc.(marker_names{lk});
        distance_here =   (markers_preproc.(marker_names{ll})-markers_preproc.(marker_names{lk}));
        delta_markers(ll,lk,:,4) =sqrt(sum((distance_here).^2,2));
        
        
       % delta_markers_reshaped = cat(2,  delta_markers_reshaped,markers_preproc.(marker_names{ll})-markers_preproc.(marker_names{lk}));
        
    end
end


if (do_video_analysis)
    %loop over different possibilities fro number ofinitial frames missing
    %to see which improves fit
    video_fps = 60;
    %video_missingframes = floor((50/60)*fps); %note this is in terms of the normal traces
    video_pc_traces_resampled = cell(1,numel(videofeatures_directory));
    inds_to_fit = cell(1,numel(videofeatures_directory));
    video_lengths = zeros(1,numel(videofeatures_directory));
    
    
    
    %% load in the video features or create them if not present
    for lk = 1:numel(videofeatures_directory)
        %% get the .times file
        %% load in the video frames
        times_files = dir(strcat(videofeatures_directory{lk},filesep,'*.times'));
        
        
        
        
        %% if using the last layer or intermediate layesr
        %             if (~do_conv_features)
        %                            videofile_here = videofeatures_files{lk};
        %          videofeatures = (load(videofile_here));
        %             videofeatures_here = squeeze(gather(videofeatures.videofeatures));
        %
        % else
        %  tag = 'convlayer*';
        % videostart = 'videofeatures';
        %% need to standardize name
        videofilenames = dir(strcat(videofeatures_directory{lk},filesep,feature_tag));
        % videofilenames = dir(strcat(videofeatures_directory{lk},filesep,videostart,'_',tag));
        fprintf('Loading video files for camera %f \n',lk)
        videofeatures_agg = [];
        load_inds = zeros(1,numel(videofilenames));
        for ll= 1:numel(videofilenames)
            % thanks Jan Simon!
            Str = videofilenames(ll).name;
            Key   = strrep(feature_tag,'*','');
            Index = strfind(Str, Key);
            load_inds(ll) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
        end
        [~,ind_to_load] = sort(load_inds,'ascend');
        
        for ll = ind_to_load
            fprintf('loading file %s \n',videofilenames(ll).name);
            test = load(strcat(videofeatures_directory{lk},filesep,videofilenames(ll).name));
            
            fieldname_here = fieldnames(test);
            [~,fieldagg] = max(size(getfield(test,fieldname_here{1})));
            videofeatures_agg = cat(fieldagg,videofeatures_agg, getfield(test,fieldname_here{1}));
        end
        %% take precautions here and above to definte the correct dimension
        [~,fieldmax] = max(size(squeeze(videofeatures_agg)));
        if (fieldmax == 1)
            videofeatures_agg = videofeatures_agg';
        end
        
        videofeatures_here = squeeze(videofeatures_agg);
        fprintf('finished loading video files \n')
        
        original_length = size(videofeatures_here,2);
        %for screw up
        frames_to_analyze = min(original_length,70000);
        
        videofeatures_here = videofeatures_here(:,1: frames_to_analyze);
        
        %get the PCs of the features
        [coeff,score,latent,tsquared,explained] = pca(squeeze(videofeatures_here)');
        video_pc_traces = score';
        %video_pc_traces = videofeatures_here;
        %resample the pcs
        num_pcs = 150;
        
        
        %% resample at the variable framerate
        %true_video_framerate = fps*(original_length+450+video_missingframes*(video_fps./fps))./size(markers_preproc.(marker_names{1}),1);
        true_video_framerate = 59.82;
        % video_pc_traces_resampled{lk} = resample(double(video_pc_traces(1:num_pcs,:))',round(fps*100),round(true_video_framerate*100));
        
        f = fopen(strcat(videofeatures_directory{lk},filesep,times_files.name));
        float1 = fread(f,[1,100000000],'uint64');
        frame_number = numel(float1);
        % 1000.*(frame_number./(float1(end)-float1(1)));
        frame_chunk = 5000;
        
        %preallocate
        video_pc_traces_resampled{lk} = zeros( floor(frames_to_analyze*245./true_video_framerate*1.1),num_pcs);
        frames_total = 0;
        frame_current = 0;
        fprintf('resampling video files \n')
        for kk = 1:frame_chunk: frames_to_analyze
            chunk_start = kk;
            chunk_end = kk+frame_chunk-1;
            frame_chunk_here = frame_chunk;
            if (chunk_end> frames_to_analyze)
                chunk_end =  frames_to_analyze;
                frame_chunk_here = (chunk_end-chunk_start+1);
            end
            
            resample_rate = 1000.*(frame_chunk_here./(float1(chunk_end)-float1(chunk_start)));
            fprintf('chunk %f resample rate %f \n',kk,resample_rate);
            mult_factor = 200;
            resampled_value = resample(double(video_pc_traces(1:num_pcs,chunk_start:chunk_end))',...
                round(fps*mult_factor),round(resample_rate*mult_factor));
            video_pc_traces_resampled{lk}( frame_current+1:(frame_current+size(resampled_value,1)),:) = resampled_value;
            frames_total = frames_total+size(resampled_value,1);
            frame_current = frame_current+size(resampled_value,1);
            %              video_pc_traces_resampled{lk} = cat(1, video_pc_traces_resampled{lk},...
            %                  resample(double(video_pc_traces(1:num_pcs,chunk_start:chunk_end))',round(fps*mult_factor),round(resample_rate*mult_factor)));
        end
        video_pc_traces_resampled{lk} = video_pc_traces_resampled{lk}(1:frame_current,:);
        fclose(f)
        
        %% subtract the missing frames from the mocap trace
        video_length = size(video_pc_traces_resampled{lk},1);
        video_lengths(lk) = video_length;
        
        %use the PCs to learn a mapping to marker positions
        figure(2300+lk)
        imagesc(squeeze(video_pc_traces_resampled{lk})')
        
    end
    
    
    %% loop over different 'microoffsets'
    missingframes_here = [1];
    %agg_mse = cell(1,numel(missingframes_here));
    do_plots_regression = 0;
    for video_missingframes = missingframes_here;
        fprintf('starting missing frames %f \n',video_missingframes);
        for lk = 1:numel(videofeatures_directory)
            
            inds_to_fit{lk} = (video_missingframes:(video_lengths(lk) +video_missingframes-1));
        end
        
        %% obtain regressors
        if (numel(video_lengths)>1)
            [min_length_val,min_ind] = min(video_lengths);
            ind_to_fit = inds_to_fit{min_ind};
            regression_input = video_pc_traces_resampled{1}(1:min_length_val,:);
            for lk = 2:numel(video_lengths)
                regression_input = cat(2, regression_input,...
                    video_pc_traces_resampled{lk}(1:min_length_val,:));
            end
            
        else
            ind_to_fit = inds_to_fit{1};
            regression_input = video_pc_traces_resampled{1};
        end
        
        %% save regression input
        
        % choose features to cluster over
        
        %do individual regressions for each relative marker position
        %Mdl = fitrsvm(video_pc_traces_resampled, delta_markers_reshaped(ind_to_fit,1),'Kernelfunction','rbf' );
        marker_scan = 1:numel(marker_names);
        num_markers = numel(marker_scan);
        num_poses = 1;
        do_regression_markers = 1;
        
        %% create better output features
        agg_features = [];
        % feature_mapping = cell(1,numel(marker_scan)*3);
        for ll = marker_scan
            agg_features = cat(1,agg_features,markers_preproc.(marker_names{ll})(:,:)');
        end
        [feature_coeff,feature_score,latent,tsquared,explained] = pca(agg_features','Centered','off');
        num_features_train = size(feature_score,2);
        %% break up training into different poses, based on ind to fit
        pose_indicies = cell(1,num_poses);
        headpos = markers_preproc.(marker_names{1})(ind_to_fit,:);
        headpos = bsxfun(@minus,headpos,median(headpos,1));
        headangle = atan2d(headpos(:,1),headpos(:,2));
        
        figure(1233)
        subplot(2,1,1)
        plot(headpos(:,1),headpos(:,2))
        subplot(2,1,2)
        plot(headangle)
        
        angle_bins = -180:360./(num_poses):180;
        
        headangle_binarized = zeros(size(headangle));
        
        for lk = 1:num_poses
            
            pose_indicies{lk} = sort(intersect(find(headangle>angle_bins(lk)),...
                find(headangle<=angle_bins(lk+1))),'ascend');
            headangle_binarized(intersect(find(headangle>angle_bins(lk)),...
                find(headangle<=angle_bins(lk+1)))) = lk;
        end
        num_outputs = size(markers_preproc.(marker_names{1}),2);
        
        if (do_regression_markers)
            Mdl = cell(num_markers,num_outputs,num_poses);
            fitinfo = cell(num_markers,num_outputs,num_poses);
            model_mse = cell(2,num_markers,num_outputs,num_poses);
        else
            Mdl = cell(num_features_train,1,num_poses);
            fitinfo = cell(num_features_train,1,num_poses);
            model_mse = cell(2,num_features_train,1,num_poses);
        end
        %don't overwrite these
        regression_input_original = regression_input;
        ind_to_fit_original = ind_to_fit;
        
        movies = cell(2,2,num_poses); %for real/regressed one for input, one for output
        
        for lk = 1:num_poses
            fprintf('for pose %f \n',lk);
            regression_input = regression_input_original;
            ind_to_fit = ind_to_fit_original;
            % restrict to a given pose
            regression_input = regression_input(pose_indicies{lk},:);
            ind_to_fit = ind_to_fit(pose_indicies{lk});
            
            input_length = size(regression_input,1);
            training_data = [1:floor(input_length *0.4),(floor(input_length*0.5)):input_length];
            held_out_data = (1+floor(input_length *0.4)):(floor(input_length*0.5));
            kfold_training_data= 1:floor(input_length);
            
            regression_output = markers_preproc.(marker_names{1})(ind_to_fit,:);
            kfold_predicted_score = zeros(numel(kfold_training_data) ,num_features_train);
            traintest_predicted_score = cell(1,2);
            
            do_kfold = 1;
            
            if (~do_regression_markers)
                %score*feature_coeff'
                for jj = 1:num_features_train
                    fprintf('training for feature %f \n',jj)
                    regression_output = feature_score(ind_to_fit,jj);
                    
                    if (do_kfold)
                        Mdl{jj,1,lk} = fitrlinear(regression_input(kfold_training_data,:), squeeze(regression_output(kfold_training_data)),...
                            'Learner','leastsquares','Solver','sparsa','GradientTolerance',1e-6,'BetaTolerance',1e-6,'Crossval','on' ...
                            );
                        model_mse(1,jj,1,lk) = kfoldLoss(Mdl{jj,1,lk});
                        % if (do_plots_regression)
                        figure(2000+jj+10000*lk)
                        %   subplot(2,num_outputs,1)
                        reduce_plot(Mdl{jj,1,lk}.kfoldPredict ,'+r')
                        hold on
                        reduce_plot(regression_output(kfold_training_data),'+k')
                        hold off
                        kfold_predicted_score(:,jj,1) = Mdl{jj,1,lk}.kfoldPredict;
                        
                    else
                        %                     Mdl{jj,1,lk} = fitrlinear(regression_input(training_data ,:), squeeze(regression_output(training_data )),...
                        %                         'Learner','leastsquares','Solver','sparsa','GradientTolerance',1e-6,'BetaTolerance',1e-6...
                        %                         );
                        %
                        Mdl{jj,1,lk} = fitrtree( regression_input(training_data ,:), squeeze(regression_output(training_data )));
                        
                        model_mse(1,jj,1,lk) = sqrt(Mdl{jj,1,lk}.loss(regression_input(training_data,:),regression_output(training_data,1),'LossFun'   ,'mse'));
                        model_mse(2,jj,1,lk) = sqrt(Mdl{jj,1,lk}.loss(regression_input(held_out_data,:),regression_output(held_out_data,1),'LossFun'   ,'mse'));
                        prediction{ll} = Mdl{jj,1,lk}.predict(regression_input(training_data,:));
                        
                        traintest_predicted_score{1} = cat(2, traintest_predicted_score{1},Mdl{jj,1,lk}.predict(regression_input(training_data,:)));
                        traintest_predicted_score{2} = cat(2, traintest_predicted_score{2},Mdl{jj,1,lk}.predict(regression_input(held_out_data,:)));
                        
                        %kFoldPredict(regression_input(training_data,:));
                        if (do_plots_regression)
                            figure(1000+jj+10000*lk)
                            subplot(2,1,1)
                            % reduce_plot(test ,'+r')
                            
                            reduce_plot(prediction{ll} ,'+r')
                            hold on
                            reduce_plot(regression_output(training_data),'+k')
                            hold off
                            
                            prediction{ll} = Mdl{jj,1,lk}.predict(regression_input(held_out_data,:));
                            
                            subplot(2,1,2)
                            reduce_plot(prediction{ll},'+r')
                            hold on
                            reduce_plot(regression_output(held_out_data),'+k')
                            hold off
                        end
                    end
                end
                
                
                reconstructed_pose_true  = feature_score*feature_coeff'; %agg_features';%
                
                for mm =1:2
                    
                    if (do_kfold)
                        reconstructed_pose = kfold_predicted_score*feature_coeff';
                    else
                        reconstructed_pose = traintest_predicted_score{mm}*feature_coeff';
                    end
                    
                    
                    if (mm ==1)
                        data_use_here = training_data;
                    else
                        data_use_here = held_out_data;
                    end
                    
                    %                     for jj = marker_scan
                    %                         markers_preproc_true.(marker_names{jj}) = markers_preproc.(marker_names{jj})(ind_to_fit(data_use_here),:);
                    %                         test = zeros(numel(data_use_here),num_outputs);
                    %                         for ll = 1:num_outputs
                    %                             test(:,ll) = Mdl{jj,ll,lk}.predict(regression_input(data_use_here,:));
                    %                         end
                    %                         markers_preproc_predicted.(marker_names{jj}) = test;
                    %                         %markers_preproc.(marker_names{jj}) = markers_preproc.(marker_names{jj})(ind_to_fit(held_out_data));
                    %                     end
                    %
                    reconstructed_markers = markers_preproc;
                    reconstructed_markers_true = markers_preproc;
                    
                    for zz =marker_scan
                        for running_ind = 1:3
                            reconstructed_markers.(marker_names{zz})(1:numel(data_use_here),(running_ind)) =reconstructed_pose(1:numel(data_use_here),3*(zz-1)+running_ind);
                            reconstructed_markers_true.(marker_names{zz})(1:numel(data_use_here),(running_ind)) =reconstructed_pose_true(ind_to_fit(data_use_here),3*(zz-1)+running_ind);
                        end
                    end
                    
                    
                    matlab_fr = 1;
                    frames_to_use = 25000;
                    trace_to_use = 1:numel(data_use_here);
                    frame_inds = trace_to_use(1:matlab_fr:min(frames_to_use,numel(trace_to_use)))';
                    save_movie = 1;
                    figure(370)
                    %initialize movies
                    movies{1,mm,lk}(1) = getframe(gcf);
                    movies{2,mm,lk}(1) = getframe(gcf);
                    
                    movies{1,mm,lk} = animate_markers(reconstructed_markers,frame_inds,marker_names(marker_scan),markercolor,links,movies{1,mm,lk},save_movie);
                    movies{2,mm,lk} = animate_markers(reconstructed_markers_true,frame_inds,marker_names(marker_scan),markercolor,links,movies{2,mm,lk},save_movie);
                    
                end
                
                
            else
                % specify_num_workers(8);
                %parpool('local',4)
                for jj =9%:numel(marker_names)
                    fprintf('for marker %f \n',jj);
                    %regression_input = cat(2,video_pc_traces_resampled(:));
                    regression_output = markers_preproc.(marker_names{jj})(ind_to_fit,:);
                    %regression_input = cat(2,regression_input, reshape(marker_position(setxor(1:size(marker_position,1),jj),ind_to_fit,:),numel(ind_to_fit),[]) );
                    prediction = cell(1,num_outputs);
                    
                    
                    for ll = 1%2:num_outputs
                        hyperopts = struct('AcquisitionFunctionName','expected-improvement-plus');
                        Lambda = logspace(-5,-1,15);
                        % [Mdl{jj,ll},fitinfo{jj,ll}] = fitrlinear(regression_input(training_data,:), squeeze(regression_output(training_data,ll)),...
                        %    'Lambda',Lambda, 'Learner','leastsquares','Solver','sparsa','Regularization','lasso','GradientTolerance',1e-6,'BetaTolerance',1e-6 ...
                        %      );
                        %                         [Mdl{jj,ll,lk},fitinfo{jj,ll,lk}] = fitrlinear(regression_input(training_data,:), squeeze(regression_output(training_data,ll)),...
                        %                             'Learner','leastsquares','Solver','sparsa','GradientTolerance',1e-6,'BetaTolerance',1e-6 ...
                        %                             );
                        fprintf('starting cross validated regression for output: %f marker %f \n',ll,jj);
                        if (do_kfold)
                            Mdl{jj,ll,lk} = fitrtree( [regression_input(kfold_training_data(1:100000) ,:)], squeeze(regression_output(kfold_training_data(1:100000),ll )),...
                                'KFold',5 );
                            
                               Mdl_nocv = fitrtree( [regression_input(kfold_training_data(1:100000) ,:)], squeeze(regression_output(kfold_training_data(1:100000),ll )),...
                                   'MaxNumSplits',20,'MinLeafSize',10,'KFold',5);
                               
                               loss( Mdl_nocv,regression_input(kfold_training_data(1:100000) ,:),squeeze(regression_output(kfold_training_data(1:100000),ll )),'LossFun','mse')
                               loss(Mdl_nocv,regression_input(kfold_training_data(100000:110000) ,:),squeeze(regression_output(kfold_training_data(100000:110000),ll )),'LossFun','mse')

                           Md
                            %,'MaxNumSplits',7 MinLeafSize, or MinParentSize
                            %Mdl =
                            %
                             Mdl_test = TreeBagger(100,regression_input(kfold_training_data([1:100000 120000:250000]) ,:),...
                                 squeeze(regression_output(kfold_training_data([1:100000 120000:250000]),ll )),'Method','regression','OOBPrediction','on'); 
                             
                             Mdl_test.oobError()
                            Mdl_test.error(regression_input(kfold_training_data([1:100000 120000:250000]) ,:),...
                                 squeeze(regression_output(kfold_training_data([1:100000 120000:250000]),ll )));
                             
                                Mdl_test.error(regression_input(kfold_training_data(100000:110000) ,:),...
                                 squeeze(regression_output(kfold_training_data(100000:110000),ll )));
                             
                              testval = Mdl_test.predict(regression_input(kfold_training_data(1:100000) ,:));
                            sqrt(mean( ( squeeze(regression_output(kfold_training_data(1:100000),ll ))-testval).^2))
                             
                            testval = Mdl_test.predict(regression_input(kfold_training_data(105000:115000) ,:));
                            sqrt(mean( ( squeeze(regression_output(kfold_training_data(105000:115000),ll ))-testval).^2))
                            %                                                   Mdl_test = fitrtree( [headangle_binarized regression_input(kfold_training_data ,:)], squeeze(regression_output(kfold_training_data )),...
                            %                                                       'KFold',5,'CategoricalPredictors',1 );
                            
                            model_mse{1,jj,ll,lk} = sqrt(Mdl{jj,ll,lk}.kfoldLoss('lossfun'   ,'mse'));
                            %model_mse(2,jj,ll,lk) = sqrt(Mdl{jj,ll,lk}.kfoldLoss('lossFun'   ,'mse'));
                            prediction{ll} = Mdl{jj,ll,lk}.kfoldPredict();
                            
                        else
                            Mdl{jj,ll,lk} = fitrtree( regression_input(training_data ,:), squeeze(regression_output(training_data )));
                            %model_mse{1,jj,ll,lk} = sqrt(Mdl{jj,ll,lk}.loss(regression_input(training_data,:),regression_output(training_data,ll),'LossFun'   ,'mse'));
                            %  model_mse{2,jj,ll,lk} = sqrt(Mdl{jj,ll,lk}.loss(regression_input(held_out_data,:),regression_output(held_out_data,ll),'LossFun'   ,'mse'));
                            prediction{ll} = Mdl{jj,ll,lk}.predict(regression_input(training_data,:));
                        end
                        %
                        %
                        % [Mdl{jj,ll}] = fitrsvm(regression_input(1:10:30000,:), squeeze(regression_output(1:10:30000,ll)),...
                        %     'KernelFunction','rbf','CrossVal' );
                        
                        
                        
                        %kFoldPredict(regression_input(training_data,:));
                        if (do_plots_regression)
                            if (~do_kfold)
                                figure(1000+jj+10000*lk+20*video_missingframes)
                                subplot(2,num_outputs,ll)
                                % reduce_plot(test ,'+r')
                                
                                reduce_plot(prediction{ll} ,'+r');
                                hold on
                                reduce_plot(regression_output(training_data,ll),'+k');
                                hold off
                                
                                
                                prediction{ll} = Mdl{jj,ll,lk}.predict(regression_input(held_out_data,:));
                                
                                subplot(2,num_outputs,num_outputs+ll)
                                reduce_plot(prediction{ll},'+r');
                                hold on
                                reduce_plot(regression_output(held_out_data,ll),'+k');
                                hold off
                            else
                                figure(1000+jj+10000*lk+20*video_missingframes)
                                subplot(1,num_outputs,ll)
                                % reduce_plot(test ,'+r')
                                
                                reduce_plot(prediction{ll} ,'+r');
                                hold on
                                reduce_plot(medfilt1(prediction{ll},5) ,'+b');
                                reduce_plot(regression_output(kfold_training_data,ll),'+k');
                                hold off
                                
                                
                                
                            end
                        end
                        
                    end
                end
            end
            %video_missingframes = missingframes_here;
            %agg_mse{find(video_missingframes == missingframes_here)} = model_mse;
            
            %% animate prediction model
            do_animation = 1;
            if do_animation
                
                mm_scan = 1:2;
                if (do_kfold)
                    mm_scan = 1;
                end
                for mm = mm_scan
                    if (~do_kfold)
                        if (mm ==1)
                            data_use_here = training_data;
                        else
                            data_use_here = held_out_data;
                        end
                    else
                        data_use_here = kfold_training_data;
                    end
                    
                    
                    markers_preproc_true = markers_preproc;
                    markers_preproc_predicted = markers_preproc;
                    
                    for jj = marker_scan
                        markers_preproc_true.(marker_names{jj}) = markers_preproc.(marker_names{jj})(ind_to_fit(data_use_here),:);
                        test = zeros(numel(data_use_here),num_outputs);
                        for ll = 1:num_outputs
                            if (~do_kfold)
                                test(:,ll) = Mdl{jj,ll,lk}.predict(regression_input(data_use_here,:));
                            else
                                test(:,ll) = Mdl{jj,ll,lk}.kfoldPredict();
                            end
                        end
                        markers_preproc_predicted.(marker_names{jj}) = test;
                        %markers_preproc.(marker_names{jj}) = markers_preproc.(marker_names{jj})(ind_to_fit(held_out_data));
                        
                        
                    end
                    %markers_preproc.(marker_names{1})
                    %movies{1,1,lk}
                    %movies{1,1,lk} = movie;
                    matlab_fr = 10;
                    trace_to_use = 1:size(markers_preproc_true.(marker_names{jj}),1);
                                        frames_to_use = 30000;%numel(trace_to_use);
                    frame_inds = trace_to_use(1:matlab_fr:min(frames_to_use,numel(trace_to_use)))';
                    save_movie = 1;
                    figure(370)
                    %initialize movies
                    movies{1,mm,lk}(1) = getframe(gcf);
                    movies{2,mm,lk}(1) = getframe(gcf);
                    
                                        movies{2,mm,lk} = animate_markers(markers_preproc_predicted,frame_inds,marker_names(marker_scan),markercolor,links,movies{2,mm,lk},save_movie);
                    movies{1,mm,lk} = animate_markers(markers_preproc_true,frame_inds,marker_names(marker_scan),markercolor,links,movies{1,mm,lk},save_movie);
                end
                
                
                
                save_movie_animation = 1;
                if (save_movie_animation)
                    save_tags{1,1} = 'true_training_2.mp4';
                    save_tags{2,1} = 'predicted_training_2.mp4';
                    save_tags{1,2} = 'true_testing.mp4';
                    save_tags{2,2} = 'predicted_testing.mp4';
                    
                    savedirectory_subcluster =strcat(savedirectory, save_tag,num2str(lk),filesep);
                    if (~exist(savedirectory_subcluster,'dir'))
                        mkdir(savedirectory_subcluster)
                    end
                    
                    for pp = 1:2
                        for mm = 1
                            if (numel(movies{pp,mm,lk}))
                                v = VideoWriter(strcat(savedirectory_subcluster,save_tags{pp,mm}),'MPEG-4');
                                open(v)
                                writeVideo(v,movies{pp,mm,lk})
                                close(v)
                            end
                        end
                    end
                end
            end
            fprintf('for pose %f model mse for missing frames  %f \n',lk,video_missingframes);
            %    model_mse(1,:,:,lk)
            %  model_mse(2,:,:,lk)
            
        end
        save_tag_model = strcat('Tree_regression_model_kfold',num2str(do_kfold),'_pcs_',num2str(num_pcs),'_poses_',num2str(num_poses),'_markers',num2str(numel(marker_scan)),'.mat');
                save_tag_mse = strcat('Tree_regression_model_kfold',num2str(do_kfold),'_pcs_',num2str(num_pcs),'_poses_',num2str(num_poses),'_markers',num2str(numel(marker_scan)),'mse.mat');
save_directory_here = strcat(savedirectory,save_tag,filesep);
        if (~exist(save_directory_here,'dir'))
            mkdir(save_directory_here);
        end
                savedirectory_model =strcat(save_directory_here, save_tag_model);
savedirectory_mse =strcat(save_directory_here, save_tag_mse);
        save(savedirectory_mse,'model_mse')
        save(savedirectory_model,'Mdl','-v7.3')
        for lk = 1:num_poses
            fprintf('for pose %f model mse for missing frames  %f \n',lk,video_missingframes);
            % model_mse(1,:,:,lk)
            %model_mse(2,:,:,lk)
        end
        
    end
end




%get aggregate feature matrix
%% simple bad frame detector
fprintf('finding bad frames \n')
tic

 bad_1 = find((delta_markers(:,:, :,4))>400);
      [~,~,i3] = ind2sub(size((delta_markers(:,:, :,4))),bad_1);
      clear bad_1
      
        bad_2 = find((delta_markers(:,:, :,3))>400);
          [~,~,i4] = ind2sub(size((delta_markers(:,:, :,4))),bad_2);
      clear bad_2
          
        bad_frames  = (cat(1,(i3),(i4)));
        clear i3 i4


% for ll = 1:numel(marker_names)
%     for lk = (ll+1):numel(marker_names)
%         bad_1 = find(delta_markers((ll),(lk), :,4)>200);
%         bad_2 = find(squeeze(delta_markers((ll),(lk), :,3))>200);
%         bad_frames  = unique(cat(1,bad_frames,bad_1,bad_2));
%     end
% end
%bad_frames = unique(bad_frames);
bad_frames = unique((bad_frames));

gap_base = zeros(1,size(delta_markers,3));
gap_base(bad_frames) = 1;
% cc = bwconncomp(gap_base);
% labelnum = cc.NumObjects;
% 
% %loop over gaps and get lengths
% %could probably use arrayfun
% 
% for kk = 1:labelnum
%         gapstart = cc.PixelIdxList{kk}(1);
%         gapend = cc.PixelIdxList{kk}(end);
% 
%         pre_gap = gapstart-1000;
%         post_gap = gapend+1000;
%         
%         if (pre_gap<1)
%             pre_gap = 1;
%         end
%         if (post_gap>size(delta_markers,3))
%             post_gap = size(delta_markers,3);
%         end
%         gap_base(pre_gap:post_gap) = 1;
%         
%     end
% end
% 
% 
% bad_frames = (bsxfun(@plus,bad_frames,-1000:1000));
gap_base_conv = conv(gap_base,ones(1,1000)./1000,'same');
%gap_base_conv(gap_base_conv>0) = 1;
bad_frames = find(gap_base_conv >0);

bad_frames = reshape(bad_frames,1,[]);
bad_frames = unique(bad_frames);
bad_frames(bad_frames<1)=1;
bad_frames(bad_frames>size(delta_markers,3)) = size(delta_markers,3);

exclude_frames_clustering = bad_frames;
fprintf('Number of frames to exclude %f \n',numel(exclude_frames_clustering));
frame_length_here = size(delta_markers,3);

clustering_ind = setxor(1:frame_length_here,intersect(1:frame_length_here,exclude_frames_clustering));
clear bad_frames;
toc

fprintf('starting clustering \n')

%% Reduce dimensionality ,compute feature spectrogram, cluster
do_clustering_analysis = 1;

    
    % determine the subclusters
    number_of_subclusters = 5;
    cluster_objects = cell(1,number_of_subclusters);
    if (rat_hands_offset)
    cluster_markersets = {[1,2,3],[6,13,11,12],[6,14,15,16],[1,2,3,4,5,6,7,8,9,10],[5,9,10,17,18],[1:numel(marker_names)],};
    else
        %head, L arm R arm, head and hips, hips, full
            cluster_markersets = {[1,2,3],[6,11,12,13],[6,14,15,16],[1,2,3,4,5,6,7,8,9,10],[5,9,10,17,18],[1:numel(marker_names)],};

    %        cluster_markersets = {[1:numel(marker_names)],[1,2,3],[3,4,5,6,7],[5,6,7,8],[1,2,3,9,10]};

    end
    if (toy_dog)
        cluster_markersets = {[1:8],[1,2,3],[3,5],[3,7,8]};
    end
    
    %cluster on a subset of frames
    frame_length_here = fps*floor(time_fraction_clustering*size(marker_velocity,2)./(fps));
    clustering_frame_length = numel(clustering_ind);
    
    
    mm=6;
    
         cluster_marker_inds = cluster_markersets{mm};
        num_subcluster_markers = numel(cluster_marker_inds);
        
%               agg_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));
% 
%         for ll = 1:num_subcluster_markers
%             for lk = (ll+1):num_subcluster_markers
%                 agg_features = cat(1,agg_features,squeeze(delta_markers(cluster_marker_inds(ll),cluster_marker_inds(lk), clustering_ind,4))',...
%                     squeeze(delta_markers(cluster_marker_inds(ll),cluster_marker_inds(lk), clustering_ind,3))');
%             end
%         end
%     
                      agg_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));

        for ll = 1:numel(marker_names)
   % marker_position(ll,:,1:3) = markers_preproc.(marker_names{ll});
   % marker_velocity(ll,2:(end),1:3) = diff(markers_preproc.(marker_names{ll}));
   % marker_velocity(ll,1,1:3) = marker_velocity(ll,2,1:3);
    agg_features = cat(1,agg_features, squeeze(marker_velocity(ll,:,4)));
        end
        
        agg_features(agg_features>6) = 6;
    %stop here
downsample = 3;
maxframes = 300000;
frames_use = 1:downsample:maxframes;

[labels , feature_mu, feature_sigma]=WaveletCluster(agg_features(:,frames_use)',1:numel(frames_use));
    


    for mm = 1:number_of_subclusters ;%number_of_clusters
        fprintf('starting cluster number %f \n',mm)
        tic
        %% get features on a subset of markers
        %cluster on a specified set of markers
        cluster_marker_inds = cluster_markersets{mm};
        num_subcluster_markers = numel(cluster_marker_inds);
        
        % choose features to cluster over
        agg_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));
        
        for ll = 1:num_subcluster_markers
            for lk = (ll+1):num_subcluster_markers
                agg_features = cat(1,agg_features,squeeze(delta_markers(cluster_marker_inds(ll),cluster_marker_inds(lk), clustering_ind,4))',...
                    squeeze(delta_markers(cluster_marker_inds(ll),cluster_marker_inds(lk), clustering_ind,3))');
            end
        end
        
        %whiten the features
        agg_features = bsxfun(@rdivide,bsxfun(@minus,agg_features,mean(agg_features,2)),std(agg_features,[],2));
        agg_features(abs(agg_features)>10) =0; 
        %could do a 'pre-pca' stage
        [coeff,score,latent,tsquared,explained] = pca(agg_features');
     %   num_pcs_1 = size(agg_features,1);
        num_pcs_2 = 20;
        num_pcs_1 = 50;
        
        pc_traces = score';
        num_pcs_1 = min(num_pcs_1,size(agg_features,1));
        
        pc_spectrograms = cell(1,num_pcs_1);
        
        if ~size(score,2)
            fprintf('UH OH ERROR IN PCA \n')
        else
        
        for k=1:num_pcs_1
          %  [~,~,time_clustering,pc_spectrograms{k}] = spectrogram(agg_features(k,:),clustering_window,clustering_overlap,fps,fps);
            [~,~,time_clustering,pc_spectrograms{k}] = spectrogram(pc_traces(k,:),clustering_window,clustering_overlap,fps,fps);
            
%         
%             parameters.numPeriods = 50;
%             parameters.samplingFreq = 245;
%             parameters.maxF = 100;
%             parameters.minF = 1;
%             parameters.omega0 = 5;
%             
%             [pc_spectrograms{k},f,scales] = findWavelets(agg_features(k,:)',20,parameters);
%             pc_spectrograms{k} = pc_spectrograms{k}';
        end
        agg_spectrograms = cell2mat(pc_spectrograms'); %second dimension is time base
        
        timebase = size(pc_spectrograms{1},2);
        
        frames_per_bin = floor(size(agg_features,2)./timebase);
        
        [coeff2,score2,latent2,tsquared2,explained2] = pca(agg_spectrograms(:,:)');
        pca_agg_spectrogram = score2';
        
        
        %% do a GMM clustering
        obj = fitgmdist((pca_agg_spectrogram(1:num_pcs_2 ,:))',num_clusters,...
            'Start', 'plus', 'Options', statset('MaxIter',1000),'RegularizationValue',0.1);
        clusterobj = cluster(obj,(pca_agg_spectrogram(1:num_pcs_2 ,:))');
        cluster_objects{mm}= clusterobj;
        toc
        
        if (do_tsne)
            mappedX = tsne(pca_agg_spectrogram(1:num_pcs_2 ,:)');
            colormap = hot(50);
            figure(31)
            hold on
            for kk =1:num_clusters
                plot(mappedX(find(clusterobj==kk),1),mappedX(find(clusterobj==kk),2),'+')
            end
        end
        end
    end
    
    
    
    
    
    
            %% visualize the result of clustering

    for mmm= 1:number_of_subclusters
                    savedirectory_subcluster =strcat(savedirectory, save_tag,num2str(mmm),filesep);
mkdir(savedirectory_subcluster)
        % figure(30)
        % plot3(squeeze(coeff(:,1)'*agg_features),squeeze(coeff(:,2)'*agg_features),squeeze(coeff(:,3)'*agg_features),'o')
        
        if (toydog)
            %     x_axis = 0:1/fps:(marker_frame_length-1)./fps;
            %
            %     figure(35)
            %     plot(x_axis,squeeze(delta_markers(3,6,clustering_ind,4)),'r')
            %     hold on
            %     plot(x_axis,squeeze(delta_markers(2,3,clustering_ind,4)),'b')
            %     area((1:frames_per_bin:frames_per_bin*numel(cluster_objects{mmm}))./fps,cluster_objects{mmm})
            %     hold off
            %
            %
            %     figure(37)
            %     plot(x_axis,squeeze(marker_velocity(2,clustering_ind,4)),'r')
            %     hold on
            %     plot(x_axis,squeeze(marker_velocity(3,clustering_ind,4)),'b')
            %     area((1:frames_per_bin:frames_per_bin*numel(cluster_objects{mmm}))./fps,0.1*cluster_objects{mmm})
            %     hold off
        end
        
        %% make spectrograms (same for all sub-marker sets)
        [~,~,times_here,spectrogram_here] = spectrogram(squeeze(delta_markers(3,6,clustering_ind,4)),clustering_window,clustering_overlap,fps,fps);
        [~,~,~,spectrogram_here_2] = spectrogram(squeeze(delta_markers(3,8,clustering_ind,4)),clustering_window,clustering_overlap,fps,fps);
        [~,~,~,spectrogram_here_3] = spectrogram(squeeze(delta_markers(5,6,clustering_ind,4)),clustering_window,clustering_overlap,fps,fps);
          [~,~,~,spectrogram_here_4] = spectrogram(squeeze(markers_preproc.(marker_names{1})(clustering_ind,3)),clustering_window,clustering_overlap,fps,fps);

%                [~,~,~,spectrogram_here_4] = spectrogram(marker_velocity(13,:,4),clustering_window,clustering_overlap,fps,fps);

        
        %  [~,~,~,spectrogram_here] = spectrogram(squeeze(delta_markers(2,3,1:frame_length_here,3)),floor(fps),floor(fps./2),fps,fps);
        time_ordering = [];
        cluster_vals = [];
        time_ordering_fulltrace = cell(1,num_clusters);
        cluster_ordering_fulltrace = cell(1,num_clusters);
      
        pose_ordering = cell(numel(marker_names),numel(marker_names));
        pose_ordering_full = cell(numel(marker_names),numel(marker_names));
        time_axis_full =  zeros(1,size(delta_markers,3));
        cluster_size =zeros(1,num_clusters);
        
            %% initialize 

        for jj = 1:numel(marker_names)
            for lk = (jj+1):num_markers
                pose_ordering_full{jj,lk} = zeros(1,size(delta_markers,3));
                pose_ordering{jj,lk}
            end
        end
        head_height = zeros(1,size(delta_markers,3));
        
        %% get positions in xyz for markers in each cluster
        fprintf('getting cluster temporal positions \n')
        lastind = zeros(num_markers,num_markers);
        lastclust = zeros(num_markers,num_markers);
        cluster_quality = zeros(3,num_clusters);
            cluster_size = zeros(1,num_clusters);
    
      
            
        for ll = 1:num_clusters
            %% ordering of the spectrogram
            fprintf('cluster %f \n',ll)
            time_ordering = cat(1,time_ordering,find(cluster_objects{mmm}==ll));
            cluster_vals = cat(1,cluster_vals,ll*ones(numel(find(cluster_objects{mmm}==ll)),1));
            cluster_size(ll) = numel(find(cluster_objects{mmm}==ll));
            
            time_ordering_fulltrace{ll} = (unique(bsxfun(@plus,round(times_here(find(cluster_objects{mmm}==ll))*clustering_window)',...
                -floor(clustering_overlap):floor(clustering_overlap))));
            
           %  cluster_ordering_fulltrace{ll} = repmat(
             
            % (bsxfun(@plus,round(times_here(find(cluster_objects{mmm}==ll))*clustering_window)',...
            %    -floor(clustering_overlap):floor(clustering_overlap)));
            
            time_ordering_fulltrace{ll}(time_ordering_fulltrace{ll}<1) = 1;
            
                        %time_ordering_fulltrace{ll} = time_ordering_fulltrace{ll}+min( clustering_ind);

            time_ordering_fulltrace{ll}(time_ordering_fulltrace{ll} > max( clustering_ind)) = max( clustering_ind);
          time_ordering_fulltrace{ll}=   clustering_ind(time_ordering_fulltrace{ll});
                %% ordering of the time series
                    for mm = find(cluster_objects{mmm}==ll)'
                        times_now = (bsxfun(@plus,round(times_here(mm)*clustering_window),(-floor(clustering_overlap):floor(clustering_overlap))));
                        times_now(times_now<1) = 1;
                      %   times_now = times_now+min( clustering_ind);
                        times_now(times_now>max( clustering_ind)) = max( clustering_ind);
                        times_now = clustering_ind(times_now);
            for jj = 1:numel(marker_names)
                for lk = (jj+1):num_markers
                    
                    
                
                      %  pose_ordering{jj,lk} = cat(1,pose_ordering{jj,lk},...
                       %     median(squeeze(delta_markers(jj,lk,times_now,4))'));
                        pose_ordering_full{jj,lk}(lastind(jj,lk)+1:lastind(jj,lk)+numel(times_now)) = (((delta_markers(jj,lk,times_now,4))));
                        head_height(lastind(jj,lk)+1:lastind(jj,lk)+numel(times_now)) = markers_preproc.(marker_names{1})(times_now,3);
                        if (jj == 1 && lk == 2)
                            time_axis_full(lastind(jj,lk)+1:lastind(jj,lk)+clustering_window) = (lastclust(jj,lk)+1)*ones(1,clustering_window);
                        end
                        lastclust(jj,lk) = lastclust(jj,lk)+1;
                        lastind(jj,lk) = lastind(jj,lk)+fps;
                        
                    end
                end
                
                
            end
            
           cluster_quality(1,ll) =  nanmean(pdist(log10(spectrogram_here(1:25,find(cluster_objects{mmm}==ll)))'));
            cluster_quality(2,ll) =  nanmean(pdist(log10(spectrogram_here_2(1:25,find(cluster_objects{mmm}==ll)))'));
             cluster_quality(3,ll) = nanmean(pdist(log10(spectrogram_here_3(1:25,find(cluster_objects{mmm}==ll)))'));
           cluster_size(1,ll) =  numel(find(cluster_objects{mmm}==ll));

        end
        
        if (toy_dog)
            xlimits = [1 500];
            xlimits = [1 numel(cluster_vals)];
            uper_freq = 75;
            figure(360+mmm)
            subplot(4,2,1)
            imagesc(log10(spectrogram_here(:,time_ordering)))
            title('head to spine')
            caxis([-8 1])
            ylim([0 50])
            xlim(xlimits)
            ylabel('Frequency [Hz]')
            
            subplot(4,2,3)
            imagesc(log10(spectrogram_here_2(:,time_ordering)))
            title('spine to foot')
            caxis([-8 0])
            ylim([0 50])
            ylabel('Frequency [Hz]')
            xlabel('time ')
            xlim(xlimits)
            
            subplot(4,2,5)
            imagesc(log10(spectrogram_here_3(:,time_ordering)))
            title('spine to tail')
            caxis([-8 1])
            ylim([0 50])
            xlim(xlimits)
            ylabel('Frequency [Hz]')
            
            subplot(4,2,7)
            area(cluster_vals')
            %xlim([1 numel(cluster_vals)])
            ylabel('Cluster Number')
            xlim(xlimits)
            title(strcat('Clustering based on markers: ',num2str(cluster_markersets{mmm})));
            
            %print('-depsc','E:\Bence\Data\MOCAP\Toydog_clusters\Clusters_20_plot.eps')
            
            %  figure(366)
            axis_ordered = 0:1./fps:(numel(pose_ordering_full{2,3})./fps-1./fps);
            subplot(4,2,2)
            plot(axis_ordered ,pose_ordering_full{2,3})
            ylim([0 100])
            
            % plot(pose_ordering{2,3})
            title('head to spine')
            %  xlim([1 numel(cluster_vals)])
            xlim(xlimits)
            ylabel('marker distance (mm)')
            
            subplot(4,2,4)
            plot( axis_ordered,pose_ordering_full{3,8})
            ylim([0 100])
            %  plot(time_axis_full,  pose_ordering_full{3,8})
            title('Spine to foot')
            % plot(pose_ordering{3,8})
            %xlim([1 numel(cluster_vals)])
            xlim(xlimits)
            ylabel('marker distance (mm)')
            ylim([0 150])
            
            subplot(4,2,6)
            plot( axis_ordered,pose_ordering_full{3,5})
            
            % plot(time_axis_full,  pose_ordering_full{3,5})
            
            % plot(pose_ordering{3,5})
            title('spine to tail')
            ylim([-150 150])
            % xlim([1 numel(cluster_vals)])
            xlim(xlimits)
            ylabel('marker distance (mm)')
            
            subplot(4,2,8)
            area(cluster_vals')
            %xlim([1 numel(cluster_vals)])
            xlim(xlimits)
            title(strcat('Clustering based on markers: ',num2str(cluster_markersets{mmm})));
            
            ylabel('Cluster Number')
            % print('-depsc',strcat(savedirectory,'toyd.eps'))
            print('-dpng',strcat(savedirectory,'toydogsubcluster',num2str(mmm),'.png'))
            
        else
           % xlimits = [1 500];
            xlimits = [1 numel(cluster_vals)];
            uper_freq = 75;
            figure(360+mmm)
            subplot(4,2,1)
            imagesc(log10(spectrogram_here(:,time_ordering)))
            title('head to spine')
            caxis([-8 1])
            ylim([0 50])
            xlim(xlimits)
            ylabel('Frequency [Hz]')
            
            subplot(4,2,3)
            imagesc(log10(spectrogram_here_2(:,time_ordering)))
            title('hip to spine')
            caxis([-8 0])
            ylim([0 50])
            ylabel('Frequency [Hz]')
            xlabel('time ')
            xlim(xlimits)
            
            subplot(4,2,5)
            imagesc(log10(spectrogram_here_4(:,time_ordering)))
            title('head height')
            caxis([-8 1])
            ylim([0 50])
            xlim(xlimits)
            ylabel('Frequency [Hz]')
            
            subplot(4,2,7)
            area(cluster_vals')
            %xlim([1 numel(cluster_vals)])
            ylabel('Cluster Number')
            xlim(xlimits)
            
            %print('-depsc','E:\Bence\Data\MOCAP\Toydog_clusters\Clusters_20_plot.eps')
            
            %  figure(366)
            axis_ordered = 0:1./fps:(numel(pose_ordering_full{3,4})./fps-1./fps);
            subplot(4,2,2)
            plot(axis_ordered ,pose_ordering_full{3,4})
            ylim([0 200])
            
            % plot(pose_ordering{2,3})
            title('head to spine')
            %  xlim([1 numel(cluster_vals)])
            xlim(xlimits)
            ylabel('marker distance (mm)')
            
            subplot(4,2,4)
            plot( axis_ordered,pose_ordering_full{3,8})
            ylim([0 200])
            %  plot(time_axis_full,  pose_ordering_full{3,8})
            
            % plot(pose_ordering{3,8})
            title('head to hip')
            %xlim([1 numel(cluster_vals)])
            xlim(xlimits)
            ylabel('marker distance (mm)')
            ylim([0 250])
            
            subplot(4,2,6)
            plot( axis_ordered,pose_ordering_full{4,5})
                       plot( axis_ordered,head_height)
 
            % plot(time_axis_full,  pose_ordering_full{3,5})
            
            % plot(pose_ordering{3,5})
            title('head height')
            ylim([-250 250])
            % xlim([1 numel(cluster_vals)])
            xlim(xlimits)
            ylabel('marker distance (mm)')
            
            subplot(4,2,8)
            area(cluster_vals')
            %xlim([1 numel(cluster_vals)])
            xlim(xlimits)
            
            ylabel('Cluster Number')
            print('-depsc',strcat(savedirectory_subcluster,'clusters_full.eps'))
            print('-dpng',strcat(savedirectory_subcluster,'clusters_full.png'))
            
        end
        
        good_clusters = find(cluster_size>=10);%intersect(find(cluster_size<25),find(cluster_size>=10)); %,
        movie_output = cell(1,numel(good_clusters));
        
        frames_use = 10000;
        movies_to_examine = good_clusters;% [1,9,17,18,19,39,41,43,95,89,50,51,59];%good_clusters(good_clusters>47);%[83,60,3,6,28,35,59,42,71,83];%[3,6,23];
        save_movie = 1;
        %save_tags = {'1_torque','9_raiseexplore','17extend','18drop1','19rise','39walknshake','41drop2','43uppershake','95walkshake','89slowwalk','50rearinvestigate','51slowshift','59turn'};
        %save_tags = {'KX83_ketamine_wobble_walk','KX60_FullTiltWobble','KX3_wobblerear','KX6_wobblewalk_headmove',...
        %    'KX28_slow_odor_search','KX35_TurnLeft','KX59_turnright','KX42wobbleup','KX71_tighttwist','KX83_fastodorsearch'};
        
        save_tags_temp = (num2str(good_clusters'));
        for kk =1:numel(good_clusters)
            save_tags{kk} = save_tags_temp(kk,:);
        end
       
        
          
        for jjj = good_clusters;%58;%36;%good_clusters(1:end)
            figure(470)
            % title(strcat('cluster ',num2str(jjj)))
      %      frame_inds = time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use,numel(time_ordering_fulltrace{jjj})));
            %          temp = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
            %  M = [];
            %   movie_output{jjj} = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
            
            trace_1 = bsxfun(@minus,squeeze(delta_markers(3,6,(time_ordering_fulltrace{jjj}),4)),...
                mean(squeeze(delta_markers(3,6,(time_ordering_fulltrace{jjj}),4))));
            
             trace_2 = bsxfun(@minus,squeeze(delta_markers(3,8,(time_ordering_fulltrace{jjj}),4)),...
                mean(squeeze(delta_markers(3,8,(time_ordering_fulltrace{jjj}),4))));
            
             trace_3 = bsxfun(@minus,squeeze(delta_markers(4,5,(time_ordering_fulltrace{jjj}),4)),...
                mean(squeeze(delta_markers(4,5,(time_ordering_fulltrace{jjj}),4))));
            
             trace_4 = bsxfun(@minus,squeeze(delta_markers(3,6,(time_ordering_fulltrace{jjj}),4)),...
                mean(squeeze(delta_markers(3,6,(time_ordering_fulltrace{jjj}),4))));
            
            [~,~,times_here,spectrogram_here_clust] = spectrogram(trace_1,...
                clustering_window,clustering_overlap,fps,fps);
            [~,~,~,spectrogram_here_2_clust] = spectrogram(trace_2,...
                clustering_window,clustering_overlap,fps,fps);
            [~,~,~,spectrogram_here_3_clust] = spectrogram(trace_3,...
                clustering_window,clustering_overlap,fps,fps);
            
                        [~,~,times_here,spectrogram_here_4_clust] = spectrogram(squeeze(marker_velocity(1,(time_ordering_fulltrace{jjj}),4)),...
                clustering_window,clustering_overlap,fps,fps);
            
            color_axis = [-8 0];
            freq_limits = [0 40];
            %xlimits = [1 500];
            % xlimits = [1 numel(cluster_vals)];
            uper_freq = 75;
            subplot(4,2,1)
            imagesc(log10(spectrogram_here_clust))
            title('head to spine')
            caxis( color_axis)
            ylim( freq_limits)
            % xlim(xlimits)
            ylabel('Frequency [Hz]')
            
            subplot(4,2,3)
            imagesc(log10(spectrogram_here_2_clust))
            title('hip to spine')
            caxis( color_axis)
            ylim( freq_limits)
           % ylabel('Frequency [Hz]')
            xlabel('time ')
            % xlim(xlimits)
            
            subplot(4,2,5)
            imagesc(log10(spectrogram_here_3_clust))
            title('inter spine')
            caxis( color_axis)
            ylim( freq_limits)
            % xlim(xlimits)
        %    ylabel('Frequency [Hz]')
          
                 subplot(4,2,7)
            imagesc(log10(spectrogram_here_4_clust))
            title('head velocity')
            caxis( color_axis)
            ylim( freq_limits)
            % xlim(xlimits)
            %ylabel('Frequency [Hz]')
          
       
            %  figure(366)
            axis_ordered = 0:1./fps:(numel(pose_ordering_full{2,3})./fps-1./fps);
            
            subplot(4,2,2)
            plot(squeeze(delta_markers(3,4,time_ordering_fulltrace{jjj},4)))
                        title('head to spine')
            ylabel('marker distance (mm)')

                        subplot(4,2,4)
            plot(squeeze(delta_markers(3,8,time_ordering_fulltrace{jjj},4)))
  title('head to hip')
       %     ylabel('marker distance (mm)')
            
            subplot(4,2,6)
           plot(squeeze(delta_markers(4,5,time_ordering_fulltrace{jjj},4)))
  title('inter spine')
          %  ylabel('marker distance (mm)')
            
               subplot(4,2,8)
           plot(245*squeeze(marker_velocity(1,time_ordering_fulltrace{jjj},4)))
  title('marker velocity')
            ylabel('head velocity (mm/s)')
            ylim([0 250])
            
                 print('-depsc',strcat(savedirectory_subcluster,'plotsfor',num2str(jjj),'.eps'))
                print('-dpng',strcat(savedirectory_subcluster,'plotsfor',num2str(jjj),'.png'))
        
%             figure(470)
%                         plot(pose_ordering_full{2,3}...
%                             (bsxfun(@plus,round(times_here(find(cluster_objects{mmm}==jjj))*clustering_window)',...
%                 -floor(clustering_overlap):floor(clustering_overlap)))')

%             ylim([0 100])
%             
%             (bsxfun(@plus,round(times_here(find(cluster_objects{mmm}==ll))*clustering_window)',...
%                 -floor(clustering_overlap):floor(clustering_overlap)))
%       
        end
        
        
         
%          
%                if (save_movie)
%             
% 
%             end
        
        for jjj = good_clusters;%58;%36;%good_clusters(1:end)
            matlab_fr = 10;
            figure(370)
           % title(strcat('cluster ',num2str(jjj)))
            frame_inds = (time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use,numel(time_ordering_fulltrace{jjj}))))';
      
           %          temp = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
           M = [];
          %  movie_output{jjj} = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
            movie_output_temp = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
           
              
                 v = VideoWriter(strcat(savedirectory_subcluster,'movie',save_tags{find(movies_to_examine==jjj)}),'MPEG-4');
                open(v)
               % writeVideo(v,movie_output{jjj})
                                writeVideo(v,movie_output_temp)

                close(v)
        end
        
        
        %             figure(370+mmm)
        %             title(num2str(jjj))
        %
        %
        %             plot3( squeeze(markers_preproc.(marker_names{1})(1,1)),...
        %                 squeeze(markers_preproc.(marker_names{1})(1,2)),...
        %                 squeeze(markers_preproc.(marker_names{1})(1,3)),'o','Color',markercolor{jj})
        %             ax = gca;
        %             axis(ax,'manual')
        %             zlim([0 250])
        %             xlim([-200 250])
        %             ylim([-200 250])
        %
        %             matlab_fr = 10;
        %
        %             frame_inds = time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use,numel(time_ordering_fulltrace{jjj})));
        %
        %             %M{jjj} = movie;
        %             %Mhere = movie;
        %             for lk = frame_inds'%1:10:10000
        %                 ind_to_plot = lk;
        %
        %
        %                 set(gca,'Nextplot','ReplaceChildren');
        %                 handles_here = cell(1,numel(marker_names));
        %                 for jj = 1:numel(marker_names)
        %                     handles_here{jj} = plot3( squeeze(markers_preproc.(marker_names{jj})(ind_to_plot,1)),...
        %                         squeeze(markers_preproc.(marker_names{jj})(ind_to_plot,2)),...
        %                         squeeze(markers_preproc.(marker_names{jj})(ind_to_plot,3)),'o','Color',markercolor{jj},'MarkerFaceColor','auto');
        %                     hold on
        %                 end
        %                 for mm = 1:numel(links)
        %                     plot3( [squeeze(markers_preproc.(marker_names{links{mm}(1)})(ind_to_plot,1)) ...
        %                         squeeze(markers_preproc.(marker_names{links{mm}(2)})(ind_to_plot,1)) ],...
        %                         [ squeeze(markers_preproc.(marker_names{links{mm}(1)})(ind_to_plot,2))...
        %                         squeeze(markers_preproc.(marker_names{links{mm}(2)})(ind_to_plot,2))],...
        %                         [squeeze(markers_preproc.(marker_names{links{mm}(1)})(ind_to_plot,3))...
        %                         squeeze(markers_preproc.(marker_names{links{mm}(2)})(ind_to_plot,3))],'Color',markercolor{links{mm}(1)});
        %                 end
        %                 %delete
        %                 if (save_movie)
        %                     M{jjj}(find(frame_inds == lk)) = getframe(gcf);
        %                 end
        %
        %                 drawnow
        %                 hold off
        %             end
        %             set(gca,'Nextplot','add');
        
        
     
            
%             for jjj = movies_to_examine;%58;%36;%good_clusters(1:end)
%                 
%                  v = VideoWriter(strcat(savedirectory_subcluster,'slow',save_tags{find(movies_to_examine==jjj)}),'MPEG-4');
%                 open(v)
%                 writeVideo(v,movie_output{jjj})
%                 close(v)
%             end
        
        plot_simul = 0;
        if (plot_simul)
           nrows = 2;
           ncols = 3;
           
         %  cluster_numbers = [2,42,51,55,59,124,73,69,65,62,183,113,112,164,158,175,15,183];
         %             cluster_numbers = [6,50,11,43,35,33,29,22,16,13,12,10];
cluster_numbers = [2,6,25,83,81,100];

     %      good_clusters(1:16);
%            cluster_names = {'tap/drop','groom hands','odor sample','up/down/up',....
%                'under and up','tap and rise','sniff and explore','left tap','lick',...
%                'tap/drop','hi sample','tap','lick and up','med sample','lever sample','low sample'};
%                     cluster_names = {'hi sample','sniff/explore','mid sniff','mid sample','low sample','head crane','reach hisample','low sniff',...
%                         'eating','sniff low','down and up','sniff and shake','sniff and explore'};    
                     cluster_names = {'CCW turn','CW turn','straight','tight turn','CCW hi tilt','CW hi tilt'};    
               
              
           %num2str(cluster_numbers');
         fighand =   figure(388)
          set(fighand,'Color','k')
         set(fighand,'Position',[100 100 1100 1100])
           for lk = 1:1000 %numel(movie_output{jjj})
           for ll = 1:nrows*ncols  %numel(cluster_numbers)
               mov_ind = cluster_numbers(ll);
               subplot_tight(nrows,ncols, ll)
               movie_size = numel(movie_output{mov_ind});
               frame_use = mod(lk,movie_size);
               if (frame_use==0)
                   frame_use = 1;
               end
               
               imshow(movie_output{mov_ind}(frame_use).cdata,movie_output{mov_ind}(frame_use).colormap)
               title(cluster_names{ll},'Color','w')
           end
             M_here(lk) =      getframe(gcf);

           end
        
v = VideoWriter(strcat(savedirectory_subcluster,'aggregate_movie_2'),'MPEG-4');
                open(v)
                writeVideo(v, M_here)
                close(v)
        clear M_here
    end
    end

    
    do_task_analysis = 1;
    if (do_task_analysis)
        
        %% get lick onsets
       speaker_on = find( abs(analog.SPEAKER)>0.3);
       speaker_on_times = zeros(1,numel(analog.SPEAKER));
       speaker_on_times(speaker_on) = 1;
     speaker_on_times = conv(speaker_on_times,ones(1,100)./100);
        speaker_on_times( speaker_on_times >0) = 1;
        
        [speaker_labels,num_tones] = bwlabel(speaker_on_times);
        
        %% get lick onsets
            lick_on = find( abs(analog.LICK_SENSOR)>0.3);
       lick_on_times = zeros(1,numel(analog.LICK_SENSOR));
       lick_on_times(lick_on) = 1;
     lick_on_times = conv(lick_on_times,ones(1,2000)./2000);
        lick_on_times( lick_on_times >0) = 1;
        
        [lick_labels,num_licks] = bwlabel(lick_on_times);
        
        
        figure(1234)
plot(analog.LICK_SENSOR)
hold on
plot(lick_on_times,'r')
plot(speaker_on_times,'g') 
plot(lever_thresholded,'c')

speaker_onsets = zeros(1,num_tones);
lick_onsets = zeros(1,num_licks);


for ll = 1:num_tones
   speaker_onsets(ll) = floor(find(speaker_labels == ll,1,'first')./20); 
      speaker_onsets(ll) = floor(find(speaker_labels == ll,1,'first')./20); 
end

for ll = 1:num_licks
       lick_onsets(ll) = floor(find(lick_labels == ll,1,'first')./20); 
      lick_onsets(ll) = floor(find(lick_labels == ll,1,'first')./20); 
end

%% get time around
   time_pre_lick = 1;
        time_post_lick = 1;
        timing_lick = -time_pre_lick*fps:time_post_lick*fps;

        good_lick_times = reshape(bsxfun(@plus,speaker_onsets',timing),1,[]);
        good_lever_times = reshape(bsxfun(@plus,speaker_onsets',timing),1,[]);

        %% get the incorrect lick onsets
        [success_licks,ind1,ind2] = intersect(lick_onsets,good_lick_times);
        lick_onsets_bad = lick_onsets;
        lick_onsets_bad(ind1) = [];
        
        %% time 
        resampled_lever_thresh = zeros(1,floor(numel(lever_thresholded)./20));
        resampled_lever_thresh(floor(find(lever_thresholded==1)./20)) = 1;
                resampled_lever_thresh(floor(find(lever_thresholded==-1)./20)) = -2;

        resampled_lever_thresh(resampled_lever_thresh>0.05) = 1;
        resampled_lever_thresh(resampled_lever_thresh<-0.05) = -1;
           resampled_lever_thresh(intersect(find(resampled_lever_thresh<0.1),...
               find(resampled_lever_thresh>-0.1))) = 0;
           figure(324)
        plot(x_axis,resampled_lever_thresh,'k')

hold on
plot(x_axis_analog,lever_thresholded,'g')
        hold off
        
        %% get the lever onsets of unsuccessful periods
           lever_onsets_bad = zeros(1,numel(lick_onsets_bad));
        for kk = 1:numel(lick_onsets_bad)
            difference_here = find(resampled_lever_thresh==1)-lick_onsets_bad(kk);
            tap_times =find(resampled_lever_thresh==1);
            [~,ind_min] = find(difference_here<0,1,'last'); %find the tap preceeding
            lever_onsets_bad(kk) = tap_times(ind_min);
        end
        
         %% get the lever onsets of successful periods
           lever_onsets_good = zeros(1,numel(speaker_onsets));
        for kk = 1:numel(speaker_onsets)
            difference_here = find(resampled_lever_thresh==1)-speaker_onsets(kk);
            tap_times =find(resampled_lever_thresh==1);
            [~,ind_min] = find(difference_here<0,1,'last'); %find the tap preceeding
            lever_onsets_good(kk) = tap_times(ind_min);
        end
        
%% get time pre post for plotting
        time_pre = 2.5;
        time_post = 1;
        timing = floor(-time_pre*fps:time_post*fps);
 %%
               time_pre_template = 1.3;
        time_post_template = 0.2;
        timing_template = floor(-time_pre_template*fps:time_post_template*fps);
        
        
reach_times = bsxfun(@plus,speaker_onsets',timing);
unsuccess_reach_times = floor(bsxfun(@plus,lick_onsets_bad',timing));
unsuccess_lever_times = floor(bsxfun(@plus,lever_onsets_bad',timing));
success_lever_times = floor(bsxfun(@plus,lever_onsets_good',timing));

unsuccess_reach_times(unsuccess_reach_times<1) = 1;
unsuccess_reach_times(unsuccess_lever_times<1) = 1;

%unsuccess_reach_times(unsuccess_reach_times<1) = 1;


lhand = markers_preproc.LHand(:,3);
lhandx = markers_preproc.LHand(:,1);
lhandy = markers_preproc.LHand(:,2);
lhand_vel = conv(marker_velocity(13,:,3),ones(1,20)./20,'same');

rhand = markers_preproc.RHand(:,3);
fhead = markers_preproc.FHead(:,3);
rspine = markers_preproc.RSpine(:,3);


lever_template = bsxfun(@plus,lever_onsets_good',timing_template);
lhandz_template = (lhand(lever_template(1,:))-mean(lhand(lever_template(1,:))))./std(lhand(lever_template(1,:)));
[xcorr_val,lagval] = xcorr(lhand,lhandz_template);

figure(1234)
plot(x_axis,0.001*xcorr_val(lagval>-1))
hold on
plot(x_axis_analog,scale_factor*analog.SPEAKER,'y')
plot(x_axis_analog,lever_thresholded,'g')

hold off
%marker_velocity(13,lever_template,4);

        fighere = figure(233)
        set(fighere,'color','w')
        subplot(1,3,1)
      plot(lhand(reach_times)')
      box off
      xlabel('Time (frames)')
      ylabel('z position (mm)')
            title('left hand')

         subplot(1,3,2)
       % plot(lhandx(success_lever_times(10:12,:))')
                plot(fhead(reach_times)')
xlabel('Time (frames)')
      ylabel('z position (mm)')
            title('front of head')
box off
            
         subplot(1,3,3)
        plot(rhand(reach_times)')
xlabel('Time (frames)')
      ylabel('z position (mm)')
      title('right hand')
      box off
                 print('-depsc',strcat(savedirectory_subcluster,'success_levertaps.eps'))
                print('-dpng',strcat(savedirectory_subcluster,'success_levertaps.png'))
      
                
        figure(234)
        subplot(1,3,1)
        plot(lhand(reach_times)')
        
        subplot(1,3,2)
          plot(rhand(reach_times)')
      
           subplot(1,3,3)
          plot(rspine(reach_times)')
          
%             subplot(1,3,3)
%           plot(fhead(reach_times)')
%       
       % hold on
        %plot(speaker_on_times,'r')


          fighere = figure(235)
        set(fighere,'color','w')
        subplot(1,3,1)
      plot(lhand(unsuccess_lever_times)')
      box off
      xlabel('Time (frames)')
      ylabel('z position (mm)')
            title('left hand')

         subplot(1,3,2)
       % plot(lhandx(success_lever_times(10:12,:))')
                plot(fhead(unsuccess_lever_times)')
xlabel('Time (frames)')
      ylabel('z position (mm)')
            title('front of head')
box off
            
         subplot(1,3,3)
        plot(rhand(unsuccess_lever_times)')
xlabel('Time (frames)')
      ylabel('z position (mm)')
      title('right hand')
      box off
                 print('-depsc',strcat(savedirectory_subcluster,'unsuccess_levertaps.eps'))
                print('-dpng',strcat(savedirectory_subcluster,'unsuccess_levertaps.png'))
      
    end
    
    
    
%% randomly remove 1% of the data for a single marker and impute repeatedly
number_of_repeats = 10;
removal_fraction = 1; %percent
buffer_length = 20;

% first select frames not within 20 frames of a true bad frame
bad_times_buffer = bsxfun(@plus,bad_times,-buffer_length:1:buffer_length);
bad_times_buffer = unique(reshape(bad_times_buffer,1,[]));
bad_times_buffer(bad_times_buffer<1) = 1;
bad_times_buffer(bad_times_buffer>marker_frame_length) = marker_frame_length;
good_times_buffered = setxor(1:marker_frame_length,bad_times_buffer);

% find times with a certain number of frames beforehand
marker_remove = 1;
random_times = randsample(good_times_buffered,floor(0.01*removal_fraction*marker_frame_length));
%build the data table for the regression
% svmtable = table();
% for j =1:numel(marker_names)
%     for lk = 1:3
%     svmtable.(strcat(marker_names{j},num2str(lk))) = markers_preproc.(marker_names{j})(:,lk);
%     end
% end
%
%
%
% Mdl = fitrsvm(svmtable(good_times_buffered,:),strcat(marker_names{marker_remove},num2str(1)));

%% plot the data
%interesting things to plot: DISTANCES spine length (3,5) shows 4 hz oscillatory
%features, also higher freq
% time series of offset to spine mid (3,4) shows 'shake' instances
% velocity of marker 4 also shows the shake

x_axis = 0:1/fps:(marker_frame_length-1)./fps;
x_axis_analog = 0:1/analog_fps:(numel(lever_thresholded)-1)./analog_fps;


scale_factor = 25;
%plot
figure(24)
%plot(x_axis,squeeze(delta_markers(2,3,:,4)),'b')
%plot(x_axis,squeeze(delta_markers(2,3,:,4)))
reduce_plot(x_axis,markers_preproc.LHand(:,3),'r')
hold on

reduce_plot(x_axis,markers_preproc.RHand(:,3),'b')
reduce_plot(x_axis,markers_preproc.FHead(:,3),'c')

%plot(x_axis,markers_preproc.(marker_names{2})(:,3),'r')

reduce_plot(x_axis_analog,scale_factor*analog.LICK_SENSOR,'g')
reduce_plot(x_axis_analog,scale_factor*lever_thresholded,'k')
reduce_plot(x_axis_analog,scale_factor*analog.SPEAKER,'y')
ylim([-50 250])
hold off

figure(30)
plot(x_axis,squeeze(delta_markers(7,8,:,4)),'r')
hold on
plot(x_axis,squeeze(delta_markers(2,3,:,4)),'b')
hold off

matrix_here = squeeze(delta_markers(3,5,:,3));
matrix_here = squeeze(marker_velocity(4,:,4));
%matrix_here = squeeze(coeff(:,8)'*agg_features);

%spectrogram with a 1s window (fps) 50% overlap (fps/2) with 1 Hz sampling
%(Nfft = fps) and the given framerate
[S,F,T,P] = spectrogram(matrix_here,fps,floor(fps./2),fps,fps);
%
% params.tapers = [2 10];
% params.Fs = 245;
% [S,t,f]= mtspecgramc(matrix_here,[1 1],params);
figure(25)
imagesc(log10(P))

figure(26)
plot(x_axis,marker_velocity(4,:,4))


figure(23)
subplot(2,1,1)
plot(squeeze(markers.R_Head))
hold on
plot(squeeze(markers.Tail))

%plot(bsxfun(@times,(resample_analog'),[1000 10^4]')')
hold off
axis tight
subplot(2,1,2)
%plot(squeeze(analog.LEVER))
plot(squeeze(markers.L_Foot))
hold on
plot(squeeze(markers.R_Foot))

%plot(squeeze(analog.SPEAKER))

axis tight

time_inds = 1:5000;
figure(24)
plot3(markers.R_Head(time_inds ,1),markers.R_Head(time_inds ,2),markers.R_Head(time_inds ,3))
hold on
plot3(markers.Tail(time_inds ,1),markers.Tail(time_inds ,2),markers.Tail(time_inds ,3),'r')
hold off

find(analog.SPEAKER==0,1,'first')