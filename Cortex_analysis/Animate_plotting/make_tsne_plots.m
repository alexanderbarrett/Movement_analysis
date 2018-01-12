function [mapped_out] = make_tsne_plots(ML_features,tsnenames)

subset = 1:100:100*7000;%size(ML_features.pose_score,1);

if ismember(tsnenames,'pose')
mappedX = tsne(cat(2,ML_features.pose_score(subset,1:6),ML_features.appearance_features_agg_score_whitened(subset,1:6)));
figure(665)
plot(mappedX(:,1),mappedX(:,2),'+')
end

mappedX_dyn_angle =  tsne(cat(2,ML_features.spectrogram_pcs_head(subset,1:15),ML_features.spectrogram_pcs_trunk(subset,1:15)));
figure(666)
plot(mappedX_dyn_angle(:,1),mappedX_dyn_angle(:,2),'+')



mappedX_dyn =  tsne(cat(2,ML_features.spectrogram_pcs_head(subset,1:5),ML_features.spectrogram_pcs_trunk(subset,1:5),...
   ML_features.spectrogram_pcs_hipL(subset,1:5),ML_features.spectrogram_pcs_hipR(subset,1:5)));
figure(666)
plot(mappedX_dyn(:,1),mappedX_dyn(:,2),'+')



mappedX_joint =  tsne(cat(2,ML_features.pose_score(subset,1:10),ML_features.appearance_features_agg_score_whitened(subset,1:6),ML_features.spectrogram_pcs_head(subset,1:10),ML_features.spectrogram_pcs_trunk(subset,1:10),...
   ML_features.spectrogram_pcs_hipL(subset,1:10),ML_features.spectrogram_pcs_hipR(subset,1:10)));

figure(667)
plot(mappedX_joint(:,1),mappedX_joint(:,2),'+')
