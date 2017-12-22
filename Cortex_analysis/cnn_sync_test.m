%% load mocap struct
%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';

mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots_videoexamples',filesep);
mkdir(savedirectory);



%% load or create struct
%createmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);
[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_videotest',mocapmasterdirectory);

[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_caff',mocapmasterdirectory);

[mocapstruct_caff] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes,0,0,[],mocapvideodirectory);

[modular_cluster_properties_social] = get_modularclustproperties(mocapstruct_social_2);


[modular_cluster_properties_social] = get_modularclustproperties(mocapstruct_social_2);

%animate_markers_aligned_fullmovie(mocapstruct_social,modular_cluster_properties_social.clustering_inds_agg{2}(1:20:end))
make_cnn_features(mocapvideodirectory)


videofeatures_agg = load('Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Recording_day7_caffeine\CameraR\636389836429298895\videofeatures_outputlayer1_R.mat');



%% If the features are saved as separate files 
%  videofilenames = dir(strcat(videofeatures_directory{lk},filesep,feature_tag));
%         % videofilenames = dir(strcat(videofeatures_directory{lk},filesep,videostart,'_',tag));
%         fprintf('Loading video files for camera %f \n',lk)
%         videofeatures_agg = [];
%         load_inds = zeros(1,numel(videofilenames));
%         for ll= 1:numel(videofilenames)
%             % thanks Jan Simon!
%             Str = videofilenames(ll).name;
%             Key   = strrep(feature_tag,'*','');
%             Index = strfind(Str, Key);
%             load_inds(ll) = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
%         end
%         [~,ind_to_load] = sort(load_inds,'ascend');
%         
%         for ll = ind_to_load
%             fprintf('loading file %s \n',videofilenames(ll).name);
%             test = load(strcat(videofeatures_directory{lk},filesep,videofilenames(ll).name));
%             
%             fieldname_here = fieldnames(test);
%             [~,fieldagg] = max(size(getfield(test,fieldname_here{1})));
%             videofeatures_agg = cat(fieldagg,videofeatures_agg, getfield(test,fieldname_here{1}));
%         end


 %% take precautions here and above to definte the correct dimension
        [~,fieldmax] = max(size(squeeze(videofeatures_agg)));
        if (fieldmax == 1)
            videofeatures_agg = videofeatures_agg';
        end
        
        videofeatures_here = squeeze(reshape(videofeatures_agg.scores_intermediate,[],size(videofeatures_agg.scores_intermediate,4)));
        fprintf('finished loading video files \n')
        
        original_length = size(videofeatures_here,2);
        %for screw up
      %  frames_to_analyze = min(original_length,70000);
        
        %videofeatures_here = videofeatures_here(:,1: frames_to_analyze);
        
        %get the PCs of the features
        [coeff,score,latent,tsquared,explained] = pca(squeeze(videofeatures_here)');
        video_pc_traces = score';
        %video_pc_traces = videofeatures_here;
        %resample the pcs
        num_pcs = 150;
        
        







  time_subset = find( mocapstruct_social_2.markers_preproc.HeadF(:,3)>160);

time_intersect = intersect(modular_cluster_properties_social.clustering_inds_agg{2},time_subset);


%% choose frames to synchronize,load videos
 v = VideoWriter(strcat(savedirectory,'sbys_movie_social'),'MPEG-4');
                
                open(v)   
                
                
                 
    filename = '\\olveczky.rc.fas.harvard.edu\Jesse\Data\Motionanalysis_captures\Vicon8\20170821\Preprocessed\social_2.mat';
    mocapstruct_social = load(filename);
    
    %load in the motion capture struct the field 'matched frames aligned'
    %defines the corresponding videoframe for each motion capture frame.
    %The video is at 50 Hz MOCAP at 300 Hz, so multiple motion capture
    %frames match onto each video frame. In actuality, there is a slight
    %mismatch that occurs every 30 minutes (540000 mocap frames) that is
    %corrected for below. You can just start with the first 540000 frames,
    %where this isn't an issue. 
    mocapstruct_caff = load('caff.mat');
    
    % change the videodirectory
    mocapstruct_caff.cameradirectory = INSERT VIDEODIRECTORY HERE
    %this selects the base frames to start visualization from
    base = 1500000;            
    %this defines a compensation for an odd shift that occurs every 540000
    %(30 minutes)
    %frames that I am still looking into. The data is synchronized enough
    %to use, or you can just start with the first 30 minutes
   offset = 780*floor(base./540000); 
   % this is a visualization of the movies. 
M= animate_markers_aligned_fullmovie_syncedvideo(mocapstruct_caff,-offset +(base:10:(base+10000)),...
    mocapstruct_caff.cameradirectory,mocapstruct_caff.matched_frames_aligned(base:10:(base+10000)));




   writeVideo(v,M)
                                    close(v)
                

%% also save an accelerated video
M_an = animate_markers_aligned_fullmovie(mocapstruct_social,modular_cluster_properties_social.clustering_inds_agg{2}((1:150:end)));

%load data
caffdata = load('E:\Dropbox\mocapdata_for_tim\caff.mat');
%this is just a compatibility issue that will be fixed
caffdata.mocap_struct.rest_frames = 1;
caffdata.mocap_struct.move_frames = 2:size(caffdata.mocap_struct.aligned_mean_position,1);

%get the indicies of the frames when all the markers are tracked
[modular_cluster_properties_caff] = get_modularclustproperties(caffdata.mocap_struct);
%look at cluster 2 - the spine, head and hips
cluster_pick =2;
%add the features to the cluster -- this splices and hi-passes all the data
%together around th
 [modular_cluster_properties_caff2] = get_clustering_features(caffdata.mocap_struct,modular_cluster_properties_caff,cluster_pick) ;

animate_markers_aligned_fullmovie(caffdata.mocap_struct,modular_cluster_properties_caff.clustering_inds_agg{cluster_pick}((1:20:end)));
%

