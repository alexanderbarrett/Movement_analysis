function [filepath,save_directory,markerinput] = function loadmocapdays(ratname,ratday)
%% Dependencies

% download btk from: https://code.google.com/archive/p/b-tk/downloads
% unzip, and addpath to matlab

% tutorial
% http://biomechanical-toolkit.github.io/docs/Wrapping/Matlab/_tutorial.html

% function list
% http://biomechanical-toolkit.github.io/docs/Wrapping/Matlab/annotated.html

%% read c3df

%filepath = 'E:\Bence\Data\MOCAP\Rat_video_datasets\20170220\Mocap\TOYDOG_videosync_long1.c3d';
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
mocapnames = {'20170725\Generated_C3D_files\'};
%mocapnames = {'20170721\Generated_C3D_files\'};

filenames_here = dir(strcat(mocapdir,mocapnames{1},'Vicon3_recording*.c3d'));
filepath_array = [];
filepath_ind = 1;
for kk = 1:numel(filenames_here)
    if (~numel(strfind(filenames_here(kk).name,'nolj')))
        filepath_array{filepath_ind} = strcat(mocapdir,mocapnames{1},filenames_here(kk).name);
        filepath_ind = filepath_ind+1;
    end
end
savedirectory = 'E:\Bence\Data\MOCAP\regression_animations\';
save_tag = 'Vicon3_18marker_amph_modularizationtest_three_';
%clustering_ind = [1:12000000];%setxor(1:frame_length_here,intersect(1:frame_length_here,exclude_frames_clustering));
