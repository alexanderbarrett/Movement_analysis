function [imagesout,oldvidobj,oldfile,fighandle] = mpfour_reader_singleframe(imagefile,framenumber,oldvidobj,oldfile,oldframe,fighandle)

%check old file to avoid new load
if ~strcmp(oldfile,imagefile)
vidObj = VideoReader(imagefile);
oldframe = 0;
else
    vidObj = oldvidobj;   
 end

% reset
oldvidobj =vidObj;
oldfile = imagefile;
%Determine the height and width of the frames.


%Create a MATLAB® movie structure array, s.
% 
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;
imagesout = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);


%% slow
%imagesout(1).cdata = vidObj.read(framenumber); 
%imagesout(1).cdata =read(vidObj,framenumber);
% if(framenumber==1)
%    fighandle= imshow(read(vidObj,framenumber));
% else
% 
% fighandle.CData =read(vidObj,framenumber);
% end
% imagesout = [];


%Read one frame at a time using readFrame until the end of the file is reached. Append data from each video frame to the structure array.

% %% works, butslowif discontinuous

% THIS IS TOO SLOW 
 frames_to_read = framenumber-oldframe;

if (frames_to_read>100)
 info = get(vidObj);
vidObj.CurrentTime = (framenumber-1)/info.FrameRate;
imagesout(1).cdata = readFrame(vidObj);
else
 for kk = 1:frames_to_read
 if (kk~= frames_to_read(end))
     readFrame(vidObj);
 else
     imagesout(1).cdata =  readFrame(vidObj);
 end
 end
end




% if hasFrame(vidObj)
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
