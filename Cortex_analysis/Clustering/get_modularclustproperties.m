
%% get frames to cluster over and exclude
function [modular_cluster_properties] = get_modularclustproperties(mocapstruct)
% determine the subclusters


%if (rat_hands_offset)
%    cluster_markersets = {[1,2,3],[6,13,11,12],[6,14,15,16],[1,2,3,4,5,6,7,8,9,10],[5,9,10,17,18],[5,9,10,17,18],[1:numel(mocapstruct.marker_names)]};
%else
    %head, L arm R arm, head and hips, hips, full
    cluster_markersets = {[1:3],[1:10,17,18],[4,7,11:13],[4,7,14:16],[6,8,9,18,19],[6,10,17,20],[6,8,9,10,17:20],...
        [1:numel(mocapstruct.markernames)]};
    % cluster_markersets = {[1,2,3],[6,11,12,13],[5,9,10,17,18],[6,14,15,16],[4,5,6,7,8],[1:numel(marker_names)]};
    subcluster_names = {'head','axial','Larm','Rarm','L leg','R leg','Both legs','Global'};
    
    % subcluster_names = {'head','Larm','Hips and knees','Rarm','Spine','Global'};
    %        cluster_markersets = {[1:numel(marker_names)],[1,2,3],[3,4,5,6,7],[5,6,7,8],[1,2,3,9,10]};
    number_of_clust = [20,50,100,20,100,100,100,200];
%end


number_of_subclusters = numel(cluster_markersets);


clustering_inds_agg = cell(1,number_of_subclusters);
badframe_inds_agg = cell(1,number_of_subclusters);
trace_length = size(mocapstruct.lever_thresholded,2);

for mmm=[1:number_of_subclusters]
    temp = zeros(1,trace_length );
    cluster_marker_inds = cluster_markersets{mmm};
    for ll = cluster_marker_inds
        temp(mocapstruct.bad_frames_agg{ll}) = 1;
    end
    %for aligned pose, also get rid of the middle spine
    temp(mocapstruct.bad_frames_agg{5}) = 1;
    
    clustering_inds_agg{mmm} = find(temp == 0);
    badframe_inds_agg{mmm} = find(temp == 1);
end


modular_cluster_properties = struct();
modular_cluster_properties.cluster_markersets = cluster_markersets;
modular_cluster_properties.subcluster_names = subcluster_names;
modular_cluster_properties.badframe_inds_agg = badframe_inds_agg;
modular_cluster_properties.clustering_inds_agg = clustering_inds_agg;

