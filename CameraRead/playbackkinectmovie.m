function [ output_args ] = playbackkinectmovie( colorfilepath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%------------------------------------------------
%Play back recordings
%------------------------------------------------
disp('Construct playback objects')
colourPlayback = VideoReader(colorfilepath);

%Set colour(c) playback parameters
cFrames = colourPlayback.NumberOfFrames;
cHeight = colourPlayback.Height;
cWidth = colourPlayback.Width;

%Preallocate movie structure
colourMov(1:cFrames)=struct('cdata', zeros(cHeight,cWidth,3,'uint8'),'colormap',[]);

disp('Reading colour frames one by one')

%read one frame at a time
for k = 1:100;%cFrames
    colourMov(k).cdata=read(colourPlayback,k);
end

disp('Sizing figure for colour playback')

%Size a figure based on the video's width and height
hf1=figure;
%set(hf1,'position',[150 150 cWidth cHeight])

disp('Playing Colour recording')

%play back the movie once at the video's frame rate
movie(hf1,colourMov,1,colourPlayback.FrameRate);


end

