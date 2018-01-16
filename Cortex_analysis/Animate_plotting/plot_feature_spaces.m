function plot_feature_spaces(ML_features,subset,color_here)

%% pose axes

figure(55)
hold on
plot3(ML_features.pose_score(subset,1),ML_features.pose_score(subset,2),ML_features.pose_score(subset,3),'+','MarkerSize',1,'Color',color_here)
title('Pose space')

figure(535)
hold on
%plot3(ML_features.appearance_features(subset,1),ML_features.appearance_features(subset,2),ML_features.appearance_features(subset,3),'+','MarkerSize',1,'Color',color_here)
% 
plot3(ML_features.appearance_features_agg_score_whitened(subset,1),...
    ML_features.appearance_features_agg_score_whitened(subset,2),...
    ML_features.appearance_features_agg_score_whitened(subset,3),'+','MarkerSize',1,'Color',color_here)
title('appearance space')


figure(539)
hold on
plot3(ML_features.dyn_score(subset,1),ML_features.dyn_score(subset,2),ML_features.dyn_score(subset,3),'+','MarkerSize',1,'Color',color_here)
%plot3(log(ML_features.dyn_score(subset,2)),log(ML_features.dyn_score(subset,3)),log(ML_features.dyn_score(subset,4)),'+','MarkerSize',1)
title('dynamics')



figure(58)
%dynamics_inds = 1:100:100*23000;
hold on
plot3(ML_features.spectrogram_pcs_trunk_angle(subset,1),ML_features.spectrogram_pcs_trunk_angle(subset,2),...
    ML_features.spectrogram_pcs_trunk_angle(subset,3),'+','MarkerSize',1,'Color',color_here)
title('trunk angle dynamics')


figure(56)
%dynamics_inds = 1:100:100*23000;
hold on
plot3(ML_features.spectrogram_pcs_trunk(subset,1),ML_features.spectrogram_pcs_trunk(subset,2),...
    ML_features.spectrogram_pcs_trunk(subset,3),'+','MarkerSize',1,'Color',color_here)
title('trunk dynamics')


figure(588)
hold on
%dynamics_inds = 1:100:100*23000;
plot3(ML_features.spectrogram_pcs_head_angle(subset,1),ML_features.spectrogram_pcs_head_angle(subset,2),...
    ML_features.spectrogram_pcs_head_angle(subset,3),'+','MarkerSize',1,'Color',color_here)
title('head angle dynamics')


figure(589)
hold on
%dynamics_inds = 1:100:100*23000;
plot3(ML_features.spectrogram_pcs_head(subset,1),ML_features.spectrogram_pcs_head(subset,2),...
    ML_features.spectrogram_pcs_head(subset,3),'+','MarkerSize',1,'Color',color_here)
title('head angle dynamics')


end