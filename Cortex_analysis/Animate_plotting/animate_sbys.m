

function M = animate_sbys(mocapstruct,mocapvideodirectory,framesubset)


%% convert video file directories to mp4
suffix = {'U','L','R'};
cameratype_ind = 2;
cameradirectory = strcat(mocapvideodirectory,filesep,'Camera',suffix{cameratype_ind},filesep);
videofolders = dir(cameradirectory);
good_directories = find(cat(1, videofolders.isdir)==1);

for mm = good_directories'
    if (numel(strfind(videofolders(mm).name,'6')))
        good_folder = videofolders(mm).name;
    end
end

imageseries_folder = strcat(cameradirectory,good_folder,filesep);


%% conver the mkv files to mp4 files
fprintf('converting mkv files \n')
camfolder = imageseries_folder;%strcat(imageseries_folder,cameradirectory);
read_mkv_file(camfolder,camfolder );



%% load .times file
%.times file
times_files = dir(strcat(camfolder,filesep,'*.times'));


f = fopen(strcat(camfolder,filesep,times_files.name));
float1 = fread(f,[1,100000000],'uint64');
frame_number = numel(float1);



%% synchronize .times file -- associate each frame with a frame
offset = float1(1);
video_frames = offset+round((0:size(mocapstruct.markers_preproc.HeadF,1)-1)*(1000/300));

maxind = min(numel(video_frames),(max(framesubset)+10000));
matched_frames =arrayfun(@(x) find(video_frames(x)-float1>0,1,'last'),1:maxind,'UniformOutput', false);

bad_frames = find(cellfun(@numel,matched_frames) == 0);
matched_frames_aligned = cat(2,zeros(1,numel(bad_frames)),cell2mat(matched_frames));



%% choose frames to synchronize,load videos
% v = VideoWriter(strcat(savedirectory,'sbys_movie_social'),'MPEG-4');
                
           %     open(v)                
M= animate_markers_aligned_fullmovie_syncedvideo(mocapstruct,framesubset,...
    camfolder,matched_frames_aligned);
 %  writeVideo(v,M)
                                  %  close(v)
                

%% also save an accelerated video
%M_an = animate_markers_aligned_fullmovie(mocapstruct_social,modular_cluster_properties_social.clustering_inds_agg{2}((1:150:end)));

