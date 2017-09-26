function imagesout = mpfour_reader(imagefile,numframes)
vidObj = VideoReader(imagefile);
%Determine the height and width of the frames.

vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
%Create a MATLAB® movie structure array, s.

imagesout = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);
%Read one frame at a time using readFrame until the end of the file is reached. Append data from each video frame to the structure array.

k = 1;
if isempty(numframes)
while hasFrame(vidObj)
    imagesout(k).cdata = readFrame(vidObj);
    k = k+1;
end
else
    for k=1:numframes
         imagesout(k).cdata = readFrame(vidObj);
    k = k+1;
    end
    
end
end