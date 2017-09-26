

fullsavepath = 'E:\Bence\Data\MOCAP\Kinect\20161026\';

%Creat a VideoReader object
common_savepath = 'vicon5_task1';
colorfilepath = fullfile(fullsavepath,strcat(common_savepath,'.mj2'));
depthfilepath = fullfile(fullsavepath,strcat(common_savepath,'_depth.mj2'));
ptcloudfilepath = fullfile(fullsavepath,strcat('ptclouds\',common_savepath,'_ptcloud_'));
metadatasavepath = fullfile(fullsavepath,strcat(common_savepath,'.mat'));
avifilepath = fullfile(fullsavepath,strcat(common_savepath,'_depth.avi'));


%playbackkinectdepth( depthfilepath,[750 925],avifilepath)
playbackkinectmovie(colorfilepath);
playbackkinectdepth(depthfilepath,[750 925],avifilepath);




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



