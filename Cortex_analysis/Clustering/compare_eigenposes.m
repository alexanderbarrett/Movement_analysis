function eigenpose = compare_eigenposes(mocapstruct,modular_cluster_properties,modno,time_subset,...
    mocapstruct2,modular_cluster_properties2,time_subset2) %clipped_index,cluster_markersets,marker_struct,plotdirectory)

%% first get features -- the poses from the clipped indicies
 agg_pose_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));  
    for ll = modular_cluster_properties.cluster_markersets{modno}   
        agg_pose_features = cat(1,agg_pose_features,mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{ll})...
            (intersect(modular_cluster_properties.clipped_index{modno},time_subset),:)');       
    end
    
    clustering_ind_2 = 1:3:size(agg_pose_features,2);
    
    %% run SVD
  [U,S,V] = svd(agg_pose_features(:,clustering_ind_2)', 'econ');
    pose_mean = mean(agg_pose_features,2);
    pose_std = std(agg_pose_features,[],2);
    
    
    %% repeat for second side
     agg_pose_features2 = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));  
    for ll = modular_cluster_properties2.cluster_markersets{modno}   
        agg_pose_features2 = cat(1,agg_pose_features2,mocapstruct2.markers_aligned_preproc.(mocapstruct2.markernames{ll})...
            (intersect(modular_cluster_properties2.clipped_index{modno},time_subset2),:)');       
    end
    
    clustering_ind_2_2 = 1:3:size(agg_pose_features2,2);
    
    %% run SVD
  [U,S,V2] = svd(agg_pose_features2(:,clustering_ind_2_2)', 'econ');
    pose_mean2 = mean(agg_pose_features2,2);
    pose_std2 = std(agg_pose_features2,[],2);
    
    h=figure;
    
%% compare over markers
    for ll = 1:min([20,size(V,2),size(V2,2)])
eigenpose = reshape(pose_mean+pose_std.*V(:,ll),3,[]);
eigenpose2 = reshape(pose_mean2+pose_std2.*V2(:,ll),3,[]);

hold on
C_1    = cell(1, numel(mocapstruct.markercolor));
C_1(:) = {'w'};

C_2    = cell(1, numel(mocapstruct2.markercolor));
C_2(:) = {'r'};
%subplot(5,4,ll)


plot_eigenpose_subset(eigenpose',modular_cluster_properties.cluster_markersets{modno},C_1,mocapstruct.links,h)
hold on
plot_eigenpose_subset(eigenpose2',modular_cluster_properties2.cluster_markersets{modno},C_2,mocapstruct2.links,h)

    end  
    fig = gcf;
fig.InvertHardcopy = 'off';
title('white is pre')
    view(90,90) 
  print('-dpng',strcat(mocapstruct.plotdirectory,'Eigenpose_overhead_comp.png')) 
   view(0,0) 
  print('-dpng',strcat(mocapstruct.plotdirectory,'Eigenpose_transverse_comp.png')) 
   view(-90,0) 
  print('-dpng',strcat(mocapstruct.plotdirectory,'Eigenpose_posterior_comp.png')) 
   view(-125,30) 
  print('-dpng',strcat(mocapstruct.plotdirectory,'Eigenpose_sideview_comp.png')) 
    
  
end