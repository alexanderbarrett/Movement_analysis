function [imagesout] = MP4Reader(imagefile, framenumber)

oldfile = 'oldfile.mp4'; % ## Not needed
oldframe = 0;            % ## Not needed

% Check old file to avoid new load
if ~strcmp(oldfile,imagefile)
    vidObj   = VideoReader(imagefile);
else
    fprintf('   --- Error Happened in frame # %d\n', framenumber);
end

vidHeight = vidObj.Height;
vidWidth  = vidObj.Width;
imagesout = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'), 'colormap',[]);

% Read one frame at a time using readFrame until the end of the file is reached.
% Append data from each video frame to the structure array.

frames_to_read = framenumber-oldframe;


% Another way
%imagesout.cdata = read(vidObj, frames_to_read); % ##

if (frames_to_read>100)
    vidObj.CurrentTime = (framenumber-1)/(vidObj.FrameRate);
    imagesout.cdata = readFrame(vidObj);
    
else
    
    for kk = 1:frames_to_read
        if (kk ~= frames_to_read(end))
            readFrame(vidObj);
        else
            imagesout.cdata = readFrame(vidObj);
        end
    end
end

end
