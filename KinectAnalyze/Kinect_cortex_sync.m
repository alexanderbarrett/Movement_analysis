
fps = 245;

%fullsavepath = 'E:\Bence\Data\Motionanalysis_captures\20170406\Vicon6_rat_kinect_long1\';
fullsavepath = 'E:\Bence\Data\Motionanalysis_captures\20170407\Toydog_kinectsync1\';

%% Creat a VideoReader object
common_savepath = 'Toydog_kinectsync';
colorfilepath = fullfile(fullsavepath,strcat(common_savepath,'.mj2'));
depthfilepath = fullfile(fullsavepath,strcat(common_savepath,'_depth.mj2'));

ptcloudfilepath = fullfile(fullsavepath,strcat('ptclouds\',common_savepath,'_ptcloud_'));
metadatasavepath = fullfile(fullsavepath,strcat(common_savepath,'.mat'));
avifilepath = fullfile(fullsavepath,strcat(common_savepath,'_depth.avi'));

%filepath = 'E:\Bence\Data\Motionanalysis_captures\20170404\Vicon6_task_long2\Vicon6_task_long2.c3d';


disp('Construct playback objects')
colourPlayback = VideoReader(colorfilepath);

%Set colour(c) playback parameters
cFrames = colourPlayback.NumberOfFrames;

depthPlayback = VideoReader(depthfilepath);

%Set depth(d) playback parameters
dFrames = depthPlayback.NumberOfFrames;

dFrames
cFrames

%% load in c3d
filepath = 'E:\Bence\Data\Motionanalysis_captures\20170407\Toydog_kinectsync1\Toydog_kinectsync1.c3d';

analog_fps = fps*20;
[markers,analog,resample_analog,lever_thresholded] = readc3d_jdm(filepath,fps,analog_fps);




%% get basic information about the dataset
marker_names = fieldnames(markers);
num_markers = numel(marker_names);
%% Data cleaning
% get frames to analyze based on set conditions
marker_frame_length = size(markers.(marker_names{1}),1);
markers_preproc = markers;

x_axis = 0:1/fps:(marker_frame_length-1)./fps;
x_axis_analog = 0:1/analog_fps:(numel(lever_thresholded)-1)./analog_fps;


scale_factor = 25;
%plot
figure(24)
%plot(x_axis,squeeze(delta_markers(2,3,:,4)),'b')
hold on
%plot(x_axis,squeeze(delta_markers(2,3,:,4)))
plot(x_axis,markers_preproc.(marker_names{2})(:,3),'r')
plot(x_axis_analog,scale_factor*analog.LICK_SENSOR,'g')
plot(x_axis_analog,scale_factor*lever_thresholded,'k')
plot(x_axis_analog,scale_factor*analog.SPEAKER,'y')
hold off

playbackkinectdepth( depthfilepath,[650 975],avifilepath,1000)
%playbackkinectmovie(colorfilepath);
framesum = framesumkinectmovie(colorfilepath,3);
%playbackkinectdepth(depthfilepath,[750 925],avifilepath);
figure(35)
%plot(framesum)
hold on
plot(x_axis_analog,scale_factor*analog.KINECTSYNC,'r')
hold off


% ptclouds = dir(strcat(ptcloudfilepath,'*'));
% 
% player = pcplayer([-0.5 0.5],[-0.5 0.5],[0.5 1],...
% 	'VerticalAxis','Z','VerticalAxisDir','down');
% 
% xlabel(player.Axes,'X (m)');
% ylabel(player.Axes,'Y (m)');
% zlabel(player.Axes,'Z (m)');
% % 
% for kk = 1:numel(ptclouds)
%    ptcloud_here = pcread(strcat(ptcloudfilepath,num2str(kk),'.ply')) ;
%       view(player,ptcloud_here);
% end



