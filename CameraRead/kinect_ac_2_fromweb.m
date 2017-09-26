%------------------------------------------------
%------------------------------------------------
%Code to record kinect colour and sensor data
%using code supplied on http://www.mathworks.co.uk/help/imaq/examples/using-the-                kinect-r-for-windows-r-from-image-acquisition-toolbox-tm.html
%and http://www.mathworks.co.uk/help/imaq/examples/logging-data-to-disk.html
%------------------------------------------------
%------------------------------------------------

imaqreset %deletes any image acquisition objects that exsist in memory and uploads         all adaptors loaded by the toolbox. As a result, image acquisition hardware is reset

%------------------------------------------------
%setting up video streams
%------------------------------------------------
disp('Setting up video streams');

%Call up dicertory containing utility functions
utilpath = fullfile(matlabroot, 'toolbox', 'imaq', 'imaqdemos', 'html', 'KinectForWindows');
addpath(utilpath);
savepath = 'C:\Users\RatControl\Documents\MOCAP';
date_val = '20160803';
fullsavepath = fullfile(savepath,date_val);

if (~exist(fullsavepath))
    mkdir(fullsavepath);
end

%Create the videoinput objects for the colour and depth streams
colourVid = videoinput('kinect', 1);
%preview(colourVid);
depthVid = videoinput('kinect', 2);

% Set the depth mode to near.
srcDepth = getselectedsource(depthVid);
srcColor = getselectedsource(colourVid);
%set(srcDepth, 'DepthMode' , 'Near');
%set(srcDepth, 'CameraElevationAngle', 0);

%set backlight compensation with centre priority
%set(srcColor, 'BacklightCompensation', 'CenterPriority');

disp('Video stream set-up complete');

%------------------------------------------------
%setting up record
%------------------------------------------------

% set the data streams to logging mode and to disk
set(colourVid, 'LoggingMode', 'Disk&Memory');
set(depthVid, 'LoggingMode', 'Disk&Memory');

%Set a video timeout property limit to 50 seconds from
%www.mathworks.com/matlabcentral/answers/103543-why-do-i-receive-the-error-getdata-timed-out-before-frames-were-available-when-using-getdata-in-im
set(colourVid, 'Timeout',50);
set(depthVid, 'Timeout',50);

%Creat a VideoReader object
common_savepath = 'vicon2_test6';
colorfilepath = fullfile(fullsavepath,strcat(common_savepath,'.mj2'));
depthfilepath = fullfile(fullsavepath,strcat(common_savepath,'_depth.mj2'));
ptcloudfilepath = fullfile(fullsavepath,strcat(common_savepath,'_ptcloud'));
metadatasavepath = fullfile(fullsavepath,strcat(common_savepath,'.mat'));


colourLogfile = VideoWriter(colorfilepath, 'Motion JPEG 2000');
depthLogfile = VideoWriter(depthfilepath, 'Archival');

%configure the video input object to use the VideoWriter object
colourVid.DiskLogger = colourLogfile;
depthVid.DiskLogger = depthLogfile;

%set the triggering mode to 'manual'
triggerconfig([colourVid depthVid], 'manual');

%set the FramePerTrigger property of the VIDEOINPUT objects to 100 to
%acquire 100 frames per trigger.
set([colourVid depthVid], 'FramesPerTrigger', 100);

disp('Video record set-up complete');

%------------------------------------------------
%Initiating the aquisition
%------------------------------------------------
disp('Starting Steams');

%Start the colour and depth device. This begins acquisition, but does not
%start logging of acquired data
start([colourVid depthVid]);

pause(5); %allow time for both streams to start

disp('Triggering')
tic
%Trigger the devices to start logging of data.
trigger([colourVid depthVid]);

%Retrieve the acquired data
[colourFrameData, colourTimeData, colourMetaData] = getdata(colourVid);
[depthFrameData, depthTimeData, depthMetaData] = getdata(depthVid);


stop([colourVid depthVid])
toc
disp('Recording Complete')



%save metadata with timestamps
disp('saving metadata')
save(metadatasavepath,'colourTimeData','colourMetaData','depthTimeData','depthMetaData');

% save pointclouds
disp('saving point clouds')
imframes = size(depthFrameData,4);
for i = 1:imframes
    ptCloud = pcfromkinect(depthVid,squeeze(depthFrameData(:,:,:,i)),squeeze(colourFrameData(:,:,:,i)));
    pcwrite(ptCloud,strcat(ptcloudfilepath,'_',num2str(i)),'PLYFormat','binary');
end

player = pcplayer([-0.5 0.5],[-0.5 0.5],[0.5 1],...
	'VerticalAxis','Z','VerticalAxisDir','down');

xlabel(player.Axes,'X (m)');
ylabel(player.Axes,'Y (m)');
zlabel(player.Axes,'Z (m)');
   view(player,ptCloud);

%playbackkinectmovie(colorfilepath);
range_here = [800 950];
playbackkinectdepth(depthfilepath,range_here);

close all;