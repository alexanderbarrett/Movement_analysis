%% annotation test

%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';
mocapmasterdirectory = '\\140.247.178.37\Jesse\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots',filesep);
mkdir(savedirectory);

%% load or create struct
%createmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);

%alternate: Vicon8_caff, Vicon8_dlslesion
[descriptor_struct_1,mocapfilearray1,mocapfilestruct1,mocapvideodirectory,mocapfiletimes1] =  get_mocap_files('Vicon8','Vicon8_caff',mocapmasterdirectory);

%% write commands here to concatenate the filearrays etc. together
[mocapstruct_concatenated] = preprocess_mocap_data( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,0,0);

%[mocapstruct_post] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_2,mocapfiletimes);

%compare the eigen clusters based on specific timepoins
cluster_here = [8];
[modular_cluster_properties] = get_modularclustproperties(mocapstruct_concatenated);

%% this gives you the frames when particular m
[modular_cluster_properties] = get_clustering_features(mocapstruct_concatenated,modular_cluster_properties,cluster_here)  ;

%% to try: SPINEF, SPINEL , HEAD(aligned), Or full pose (modular set 2, inc. the hips and knees and offset)

%visualize
animate_markers_aligned_fullmovie(mocapstruct_concatenated,modular_cluster_properties.clustering_inds_agg{2}(1:20:end))

% start te gui
MCC;


