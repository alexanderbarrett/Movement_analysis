function [imagesout,oldvidobj,oldfile] = mpfour_reader_singleframe(imagefile,framenumber,oldvidobj,oldfile,oldframe)

%check old file to avoid new load
if ~strcmp(oldfile,imagefile)
vidObj = VideoReader(imagefile);
else
    vidObj = oldvidobj;
 end

% reset
oldvidobj =vidObj;
oldfile = imagefile;
%Determine the height and width of the frames.

vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
%Create a MATLAB� movie structure array, s.

imagesout = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);

%Read one frame at a time using readFrame until the end of the file is reached. Append data from each video frame to the structure array.


 
 frames_to_read = framenumber-oldframe;
 for kk = 1:frames_to_read
 if (kk~= frames_to_read(end))
     readFrame(vidObj);
 else
     imagesout(1).cdata =  readFrame(vidObj);
 end
 end
 
% THIS IS TOO SLOW vidObj.CurrentTime = (framenumber-1)/info.FrameRate;
 % info = get(vidObj);

% if hasFrame(vidObj)
%imagesout(1).cdata = readFrame(vidObj);
%end



%% old version
% 
% if isempty(numframes)
% while hasFrame(vidObj)
%     imagesout(k).cdata = readFrame(vidObj);
%     k = k+1;
% end
% else
%     for k=1:numframes
%          imagesout(k).cdata = readFrame(vidObj);
%     k = k+1;
%     end
%     
% end


end