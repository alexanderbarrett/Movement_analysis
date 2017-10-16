
%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';
mocapmasterdirectory = '\\140.247.178.37\Jesse\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots',filesep);
mkdir(savedirectory);

%% load or create struct
%createmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);

[descriptor_struct_1,mocapfilearray,mocapfilestruct] =  get_mocap_files('Vicon8','Vicon8_prelesion',mocapmasterdirectory);
[mocapstruct_pre] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1);

[descriptor_struct_2,mocapfilearray,mocapfilestruct] =  get_mocap_files('Vicon8','Vicon8_dlslesion_late',mocapmasterdirectory);
[mocapstruct_post] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_2);

plotfractionmissing(mocapstruct_pre)

compare_plot_marker_characteristics_timerange(mocapstruct_pre,mocapstruct_pre.move_frames,mocapstruct_post,mocapstruct_post.move_frames)

%compare_comp_fraction_moving(mocapstruct)
%compare_plot_marker_characteristics(mocapstruct)

%compare the eigen clusters based on specific timepoins
cluster_here = [2];
[modular_cluster_properties] = get_modularclustproperties(mocapstruct_pre);
[modular_cluster_properties2] = get_modularclustproperties(mocapstruct_post);

%% the clipped index is the index of clipped frames in the full trace
[modular_cluster_properties] = get_clustering_features(mocapstruct_pre,modular_cluster_properties,cluster_here)  ;
[modular_cluster_properties2] = get_clustering_features(mocapstruct_post,modular_cluster_properties2,cluster_here)  ;

    %% first compute the eigenposes to look for postural differences
      %% look at the eigenposes of the data
        fprintf('computing eigenposes \n')
  time_subset = find(( mocapstruct.markers_preproc.SpineF(:,3)-mocapstruct.markers_preproc.SpineL(:,3))>60);
    time_subset = find(( markers_clipped.SpineL(:,2)>20));
  
 eigenpose = compare_eigenposes(mocapstruct_pre,modular_cluster_properties,cluster_here,time_subset,...
     mocapstruct_post,modular_cluster_properties2,time_subset2);