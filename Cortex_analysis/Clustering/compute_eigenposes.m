function eigenpose = compute_eigenposes(mocapstruct,modular_cluster_properties,modno) %clipped_index,cluster_markersets,marker_struct,plotdirectory)
%% first get features -- the poses from the clipped indicies
 agg_pose_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));  
    for ll = modular_cluster_properties.cluster_markersets{modno}   
        agg_pose_features = cat(1,agg_pose_features,mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{ll})(modular_cluster_properties.clipped_index{modno},:)');       
    end
    
    clustering_ind_2 = 1:3:size(agg_pose_features,2);
    
    %% run SVD
  [U,S,V] = svd(agg_pose_features(:,clustering_ind_2)', 'econ');
    pose_mean = mean(agg_pose_features,2);
    pose_std = std(agg_pose_features,[],2);
    
    h=figure;

    for ll = 1:20
eigenpose = reshape(pose_mean+pose_std.*V(:,ll),3,[]);

hold on

%subplot(5,4,ll)
plot_eigenpose_subset(eigenpose',modular_cluster_properties.cluster_markersets{modno},mocapstruct.markercolor,mocapstruct.links,h)
hold on
    end  
    fig = gcf;
fig.InvertHardcopy = 'off';
    view(90,90) 
  print('-dpng',strcat(mocapstruct.plotdirectory,'Eigenpose_overhead.png')) 
   view(0,0) 
  print('-dpng',strcat(mocapstruct.plotdirectory,'Eigenpose_transverse.png')) 
   view(-90,0) 
  print('-dpng',strcat(mocapstruct.plotdirectory,'Eigenpose_posterior.png')) 
   view(-125,30) 
  print('-dpng',strcat(mocapstruct.plotdirectory,'Eigenpose_sideview.png')) 
    
  
end