function visualize_mlfeat_eigenposes(ML_features,mocapstruct,ep_no,h,markersetind)


C_1    = cell(1, numel(mocapstruct.markercolor));
C_1(:) = {'w'};
C_2    = cell(1, numel(mocapstruct.markercolor));
C_2(:) = {'r'};
C_3    = cell(1, numel(mocapstruct.markercolor));
C_3(:) = {'b'};

plot_eigenpose_subset(reshape(ML_features.pose_mean+std(ML_features.pose_score(:,ep_no))*ML_features.pose_coeffs(:,ep_no),3,[])',mocapstruct.modular_cluster_properties.cluster_markersets{markersetind},C_2,mocapstruct.links,h)
hold on
plot_eigenpose_subset(reshape(ML_features.pose_mean,3,[])',mocapstruct.modular_cluster_properties.cluster_markersets{markersetind},C_1,mocapstruct.links,h)
plot_eigenpose_subset(reshape(ML_features.pose_mean-std(ML_features.pose_score(:,ep_no))*ML_features.pose_coeffs(:,ep_no),3,[])',mocapstruct.modular_cluster_properties.cluster_markersets{markersetind},C_3,mocapstruct.links,h)
hold off
view(90,90) 

end
