function mocapstruct= assign_modular_annotation_properties(mocapstruct,cluster_here)


[modular_cluster_properties] = get_modularclustproperties(mocapstruct);
 [modular_cluster_properties] = get_clustering_features(mocapstruct,modular_cluster_properties,cluster_here)  ;
mocapstruct.modular_cluster_properties = modular_cluster_properties;

    annotated_pose_struct = auto_annotate_pose(mocapstruct.markers_aligned_preproc,mocapstruct);
mocapstruct.annotated_pose_struct  = annotated_pose_struct ;


end