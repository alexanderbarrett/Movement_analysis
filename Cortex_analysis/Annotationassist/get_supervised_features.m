function ML_features = get_supervised_features(mocapstruct,framelist,markersetind)


%% local pose/dynamics features
% load saved eigenposture feature coefficients and dynamics coefficients
eigenposture_save_folder = 'Y:\Jesse\Data\Annotation\Eigenposture_coefficients';
decimation_factor = 3;
%parameters for hipass clip and clustering
params.fps = 300;

opts.fps =300./1;
opts.clustering_window = opts.fps./2;
opts.clustering_overlap = opts.fps./4;
    num_eigenpcs = 10;%size(COEFF,1);
num_dynamics_pcs = 3;

pose_coeff_file = strcat(eigenposture_save_folder,filesep,'pose_coeff.mat');
dyn_coeff_file = strcat(eigenposture_save_folder,filesep,'dynamics_coeff.mat');
overwrite_coeff = 0;

framelist_true = intersect(mocapstruct.modular_cluster_properties.clipped_index{markersetind},framelist);

ML_features.framelist_true = framelist_true;

  agg_pose_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));  
    for ll = mocapstruct.modular_cluster_properties.cluster_markersets{markersetind}   
        agg_pose_features = cat(1,agg_pose_features,mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{ll})...
            (framelist_true,:)');       
    end
    clustering_ind_2 = 1:decimation_factor:size(agg_pose_features,2);

      pose_mean = mean(agg_pose_features(:, clustering_ind_2),2);
%center the pose
centered_pose = bsxfun(@minus,agg_pose_features,pose_mean);
    
% if doesnt exist compute eigenposture feature coefficients and dynamics
% coefficients
if (~exist(pose_coeff_file,'file') || overwrite_coeff)
       [COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(agg_pose_features');
    save(pose_coeff_file,'COEFF');
else
    load(pose_coeff_file);
end
%compute score
pose_score = centered_pose'*COEFF;
pose_score_whitened = pose_score;%bsxfun(@rdivide,pose_score,std(pose_score,[],1));
ML_features.pose_score = pose_score_whitened(:,1:num_eigenpcs);

%figure(44)
%plot3(pose_score(1:20:end,1),pose_score(1:20:end,2),pose_score(1:20:end,3),'+')

%whiten the eigenfeatures, if not already

%% compute the spectrogram of the whitened eigenfeatures
%get the clipped time trace of the score

   feature_mat_dyn = [];


     clipped_pre =  hipass_clip_fragments(mocapstruct.markers_aligned_preproc,framelist_true,params);
     
for mm = mocapstruct.modular_cluster_properties.cluster_markersets{markersetind}  
     feature_mat_dyn = cat(2,feature_mat_dyn,clipped_pre.(mocapstruct.markernames{mm}));
end
feature_mat_dyn_centered = bsxfun(@minus,feature_mat_dyn,mean(feature_mat_dyn,1));
         feature_mat_dyn_posespace = feature_mat_dyn_centered*COEFF;
  feature_mat_dyn_posespace_whitened = bsxfun(@rdivide,feature_mat_dyn_posespace,std(feature_mat_dyn_posespace,[],1));

    
[dyadic_spectrograms,fr,timespect] = get_dyadic_spectrogram(feature_mat_dyn_posespace(:,1:num_eigenpcs)',opts);
dyadic_spectrograms_reshaped = reshape(dyadic_spectrograms,[],num_eigenpcs,size(dyadic_spectrograms,2));
figure(45)
imagesc(squeeze(dyadic_spectrograms_reshaped(:,3,:)))
eigenpose_dynamics_coeffs = zeros(num_eigenpcs,numel(fr),num_dynamics_pcs);

% whiten each spectrogram vector
%dyadic_spectrograms_1_w = bsxfun(@rdivide,bsxfun(@minus,dyadic_spectrograms_1,mean(dyadic_spectrograms_1,2)),std(dyadic_spectrograms_1,[],2));
%get the principal components of the spectrogram

%loop over all spectrograms
if (~exist(dyn_coeff_file,'file') || overwrite_coeff)
    for ll = 1:size(dyadic_spectrograms_reshaped,2)
       [COEFF_dyn, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(squeeze(dyadic_spectrograms_reshaped(:,ll,:))');
       eigenpose_dynamics_coeffs(ll,:,:) = COEFF_dyn(:,1:num_dynamics_pcs);
    end
    save(dyn_coeff_file,'eigenpose_dynamics_coeffs');
else
load(dyn_coeff_file)
end



% compute eigenposture loadings
dyadic_spectrograms_score = [];
for ll = 1:size(dyadic_spectrograms_reshaped,2)
    dyadic_spectrograms_white = bsxfun(@minus,squeeze(dyadic_spectrograms_reshaped(:,ll,:)),mean(squeeze(dyadic_spectrograms_reshaped(:,ll,:)),2));
        dyadic_spectrograms_score_temp = dyadic_spectrograms_white'*squeeze( eigenpose_dynamics_coeffs(ll,:,:));
%dyadic_spectrograms_score_temp = bsxfun(@rdivide,dyadic_spectrograms_score_temp,std(dyadic_spectrograms_score_temp,[],2));
%plot(dyadic_spectrograms_score_temp)
    dyadic_spectrograms_score = cat(2,dyadic_spectrograms_score,dyadic_spectrograms_score_temp);
    end
ML_features.dynamics_score = dyadic_spectrograms_score;
replication_factor = floor(size(feature_mat_dyn,1)./size(dyadic_spectrograms_score,1));
ML_features.dynamics_replication_ind = reshape(repmat(1:size(dyadic_spectrograms_score,1),replication_factor,1),1,[]);
%get the frames of each spectrogrambin
%timespect_full = unique(bsxfun(@plus,timespect',-opts.clustering_window:opts.clustering_window));
%timespect_full(timespect_full<1) = 1;

%figure(45)
%plot3(dyadic_spectrograms_score(1:5:end,1),dyadic_spectrograms_score(1:5:end,4),dyadic_spectrograms_score(1:5:end,7),'+')
% compute eigenposture dynamics

%replicate array to cover the space

%% hand designed pose features
  %high rear
  ML_features.high_rear = mocapstruct.markers_aligned_preproc.SpineF(framelist_true,3)-mocapstruct.markers_aligned_preproc.SpineL(framelist_true,3);
    ML_features.very_high_rear = mocapstruct.markers_aligned_preproc.SpineF(framelist_true,3)-mocapstruct.markers_aligned_preproc.SpineL(framelist_true,3);

  %low rear -- shortens stance more and more
   ML_features.low_rear = mocapstruct.markers_aligned_preproc.HeadB(framelist_true,3)-mocapstruct.markers_aligned_preproc.SpineF(framelist_true,3);
  
  %l/r groom 
  ML_features.RGroom = mocapstruct.markers_aligned_preproc.SpineF(framelist_true,2)-mocapstruct.markers_aligned_preproc.SpineL(framelist_true,2);
    ML_features.LGroom =mocapstruct.markers_aligned_preproc.SpineF(framelist_true,2)-mocapstruct.markers_aligned_preproc.SpineL(framelist_true,2);


% get standard deviation over windows


%% morphology features -- inter marker distances
appearance_pairs = {{'HeadB','SpineF'},{'SpineF','SpineM'},{'SpineM','SpineL'},{'SpineM','Offset1'},{'Offset2','Offset1'},...
    {'Offset2','SpineL'},{'Offset1','SpineF'}};

appearance_features = zeros(size(mocapstruct.markers_preproc.HeadF,1),numel(appearance_pairs));
for ll = 1:numel(appearance_pairs)
   appearance_features(:,ll) = vectornorm( mocapstruct.markers_preproc.(appearance_pairs{ll}{1}),mocapstruct.markers_preproc.(appearance_pairs{ll}{2}),2);
end

ML_features.appearance_features = appearance_features(framelist_true,:);

%% appearance features -- joint angles
%get_segment_jointangles


%% trajectory features
% get average trunk displacement, velocity, acceleration, over different
% time windows: 10,30,100,300 frame windows
% also head head vel
difforders = [10,33];
params.difforder = 10;
params.medfiltorder = 10;
params.gaussorder = 2.5;
ML_features.trunk_vel =zeros(numel(difforders),numel(framelist_true));
ML_features.head_vel =zeros(numel(difforders),numel(framelist_true));

for ll = 1:numel(difforders)
    params.difforder = difforders(ll);
    temp1 = get_markersubgroup_velocity(mocapstruct.markers_preproc,[4:8],params);
        temp2 = get_markersubgroup_velocity(mocapstruct.markers_aligned_preproc,[1:3],params);

ML_features.trunk_vel(ll,:) = temp1(framelist_true);
ML_features.head_vel(ll,:) = temp2(framelist_true);
end

%movstd()

% figure(55)
% %plot(ML_features.trunk_vel(1,:)')
% plot(xcorr(ML_features.head_vel(2,:),2000,'coeff'))
% hold on
% plot(xcorr(ML_features.trunk_vel(2,:),2000,'coeff'),'r')
% plot(xcorr(ML_features.head_vel(2,:),ML_features.trunk_vel(2,:),2000,'coeff'),'g')
% hold off

%% window features
%get individual eigenpose velocity and standard deviation, do rolling
ML_features.pose_window_sd = zeros(numel(difforders),size(feature_mat_dyn_posespace_whitened,1),num_eigenpcs);
ML_features.pose_velocity = zeros(numel(difforders),size(feature_mat_dyn_posespace_whitened,1),num_eigenpcs);
for mm = 1:num_eigenpcs
    for ll = 1:numel(difforders)
        params.difforder = difforders(ll);
      ML_features.pose_window_sd(ll,:,mm) = movstd(feature_mat_dyn_posespace_whitened(:,mm),difforders(ll));
        ML_features.pose_velocity(ll,:,mm) = get_vector_velocity(feature_mat_dyn_posespace_whitened(:,mm),params);
       
    end
end

%% get joint angle features
saggital_inds =  [1,3];
saggital_pairs = {{'HeadB','SpineF'},{'SpineF','SpineM'} ,...
    {'SpineL','SpineM'},{'SpineF','SpineM'} };

% get appearance feature
%ML_features.eigenpose_velocity
%ML_features.eigenpose_sd


%get_filtered_derivative(pose_score_whitened(:,1),numframes);
%get_window_sd(pose_score_whitened(:,1),numframes);

% get marker relative window features
%get_markersubgroup_velocity(mocapstruct.markers_preproc,[4:8],params);
%get_filtered_derivative(pose_score_whitened(:,1),numframes);
%get_window_sd(pose_score_whitened(:,1),numframes);



end