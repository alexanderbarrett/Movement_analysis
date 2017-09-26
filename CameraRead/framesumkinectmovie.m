function [ frame_sum ] = framesumkinectmovie( colorfilepath,color_channel )
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

frame_sum = zeros(1,cFrames);
%read one frame at a time
for k = 1:cFrames
    colourMov(k).cdata=read(colourPlayback,k);
    frame_sum(k) = sum(sum(colourMov(k).cdata(:,:,color_channel)));
end


end

