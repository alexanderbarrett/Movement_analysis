function [modular_cluster_properties] = get_clustering_features(mocapstruct,modular_cluster_properties,modno) 
%% clip the makrers
    fprintf('clipping and hipass filtering marker positions \n')
    params_clip.fps = mocapstruct.fps;
    [markers_clipped,clipped_index] = hipass_clip(mocapstruct.markers_aligned_preproc,modular_cluster_properties.badframe_inds_agg{modno},params_clip);
    
    agg_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));
    
    for ll = modular_cluster_properties.cluster_markersets{modno}      
        agg_features = cat(1,agg_features,markers_clipped.(mocapstruct.markernames{ll})');
    end
    
    
    
    %% change labels for features
    feature_labels = mocapstruct.markernames;
    feature_labels = feature_labels((modular_cluster_properties.cluster_markersets{modno}));
    
    if (~isfield('modular_cluster_properties',clipped_index))
    modular_cluster_properties.agg_features = cell(1,numel(modular_cluster_properties.badframe_inds_agg));
    modular_cluster_properties.clipped_index = cell(1,numel(modular_cluster_properties.badframe_inds_agg));
        modular_cluster_properties.feature_labels = cell(1,numel(modular_cluster_properties.badframe_inds_agg));
    end
    modular_cluster_properties.agg_features{modno} = agg_features;
    modular_cluster_properties.clipped_index{modno} = clipped_index;
    modular_cluster_properties.feature_labels{modno} = feature_labels;
    
   % feature_labels_agg{mmm} = feature_labels;
    
    %% downsample for clustering
  
    
end