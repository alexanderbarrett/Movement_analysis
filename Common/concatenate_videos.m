


function outputvideo = concatenate_videos(movie_1,movie_2,concat_name)
v1 = VideoReader(movie_1);
v2 = VideoReader(movie_2);

frames1=v1.Numberofframes;
frames2=v2.Numberofframes;

min_frames = min(frames1,frames2);
dimension_combine =2;

outputvideo = VideoWriter(concat_name,'MPEG-4');
open(outputvideo);

for ll = 1:min_frames
    frame_temp = cat(dimension_combine,read(v1,ll),read(v2,ll));
    writeVideo(outputvideo,frame_temp)
end

close(outputvideo)
end
