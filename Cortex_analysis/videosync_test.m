%% load mocap struct
%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';
mocapmasterdirectory = '\\140.247.178.37\Jesse\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots',filesep);
mkdir(savedirectory);

%% load or create struct
%createmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);
[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory] =  get_mocap_files('Vicon8','Vicon8_videotest',mocapmasterdirectory);
[mocapstruct_social] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1);

[modular_cluster_properties_social] = get_modularclustproperties(mocapstruct_social);
[modular_cluster_properties_caff] = get_modularclustproperties(mocapstruct_social);

%animate_markers_aligned_fullmovie(mocapstruct_social,modular_cluster_properties_social.clustering_inds_agg{2}(1:20:end))


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
video_frames = offset+round((0:size(mocapstruct_social.markers_preproc.HeadF,1)-1)*(1000/300));

matched_frames =arrayfun(@(x) find(video_frames(x)-float1>0,1,'last'),1:1000000,'UniformOutput', false);

bad_frames = find(cellfun(@numel,matched_frames) == 0);
matched_frames_aligned = cat(2,zeros(1,numel(bad_frames)),cell2mat(matched_frames));


%% choose frames to synchronize,load videos
animate_markers_aligned_fullmovie_syncedvideo(mocapstruct_social,modular_cluster_properties_social.clustering_inds_agg{2}(600000:20:700000),...
    camfolder,matched_frames_aligned)



%% play movies side by side