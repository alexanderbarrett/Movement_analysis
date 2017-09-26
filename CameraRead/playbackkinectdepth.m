function [ output_args ] = playbackkinectdepth( depthfilepath,range_here,avifilepath,numframes )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

depthPlayback = VideoReader(depthfilepath);

%Set depth(d) playback parameters
dFrames = depthPlayback.NumberOfFrames;
dHeight = depthPlayback.Height;
dWidth = depthPlayback.Width;

%Preallocate movie structure
depthMov(1:dFrames)=struct('cdata', zeros(dHeight,dWidth,3,'uint8'),'colormap',jet(256));

disp('Reading depth frames one by one')

%read one frame at a time
maxDistFromCamera = 1600; % 1600 millimeters
VideoWriter(avifilepath)

frames_here = dFrames;
if numel(numframes)
frames_here = numframes;
end

for k = 1:frames_here
    % Depth frames are int16.
    depthFrame = read(depthPlayback,k);
   % max(depthFrame)
    depthFrame_scaled = depthFrame;
    depthFrame_scaled(depthFrame_scaled<range_here(1)) = range_here(1);
%    range_here(2)
    depthFrame_scaled(depthFrame_scaled>range_here(2)) = range_here(2);
    depthFrame_scaled = 255.0*(single(depthFrame_scaled-range_here(1))./...
        (range_here(2)-range_here(1)));
 %   max(depthFrame_scaled)
    % We'll rescale the image from [0,maxDistFromCamera] to [0,255]
    %depthFrame = 255.0*single(depthFrame)/maxDistFromCamera;
    % And then recast it to uint8 for display.
    depthMov(k).cdata=uint8(depthFrame_scaled);
    if (mod(k,100)==0)
        fprintf('on frame %f \n',k);
    end
end

disp('Sizing figure for depth playback')

%Size a figure based on the video's width and height
hf2=figure;
set(hf2,'position',[150 150 dWidth dHeight])

disp('Playing Depth recording')

%play back the movie once at the video's frame rate
playmovie = 0;
if (playmovie)
movie(hf2,depthMov,1,depthPlayback.FrameRate);
end


%% save in chunks to circumvent file size limit
chunksize = 2000;
num_chunks = ceil(frames_here./chunksize);
fprintf('number of chunks %f filepath %s',num_chunks,avifilepath)
for kk=1:num_chunks
avipathhere = strrep(avifilepath,'.avi',strcat('_',num2str(kk),'.avi'));
movie2avi(depthMov((kk-1)*chunksize+1:min(kk*chunksize,frames_here)),avipathhere,'compression','None');
end
colormap(jet)
end

