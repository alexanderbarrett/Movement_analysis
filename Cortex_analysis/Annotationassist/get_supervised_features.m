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
jadyn_coeff_file = strcat(eigenposture_save_folder,filesep,'ja_coeff.mat');
appearance_coeff_file = strcat(eigenposture_save_folder,filesep,'appearance_coeff.mat');
dyn_coeff_file_markers = strcat(eigenposture_save_folder,filesep,'dynamics_coeff_2.mat');

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
ML_features.pose_coeffs = COEFF(:,1:num_eigenpcs);
ML_features.pose_mean = pose_mean;

for ll = 1:3
figure(399)
h=subplot(1,3,ll)
visualize_mlfeat_eigenposes(ML_features,mocapstruct,ll,h,markersetind)
end

%open up memory
clear centered_pose

% 
% figure(44)
% plot3(pose_score(1:20:end,1),pose_score(1:20:end,2),pose_score(1:20:end,3),'+')

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
ML_features.eigenpose_dynamics_coeffs = eigenpose_dynamics_coeffs;


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

%open up memory
clear dyadic_spectrograms_score dyadic_spectrograms_reshaped  dyadic_spectrograms

%get the frames of each spectrogrambin
%timespect_full = unique(bsxfun(@plus,timespect',-opts.clustering_window:opts.clustering_window));
%timespect_full(timespect_full<1) = 1;
% 
% figure(45)
% plot3(dyadic_spectrograms_score(1:5:end,1),dyadic_spectrograms_score(1:5:end,2),dyadic_spectrograms_score(1:5:end,3),'o','MarkerSize',2)
%  %compute eigenposture dynamics

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
markers_to_loop = mocapstruct.modular_cluster_properties.cluster_markersets{markersetind}(1:10);
appearance_features_agg = zeros(size(mocapstruct.markers_preproc.HeadF(framelist_true,:),1),numel(markers_to_loop )*(numel(markers_to_loop )-1)./2);
marker_fn =  fieldnames(mocapstruct.markers_preproc);


%% Get the PCS of the relative distances
ind_use = 1;
for ll = 1:numel(markers_to_loop)
for jj = ll+1:numel(markers_to_loop)
   
    appearance_features_agg(:,ind_use) = ...
        vectornorm( mocapstruct.markers_aligned_preproc.(marker_fn{ll})(framelist_true,:),mocapstruct.markers_aligned_preproc.(marker_fn{jj})(framelist_true,:),2);
    ind_use = ind_use +1;
end
end
appearance_features_agg = appearance_features_agg';
appearance_mean = mean(appearance_features_agg(:,:),2);
%center the pose
centered_appearance = bsxfun(@minus,appearance_features_agg,appearance_mean);

if (~exist(appearance_coeff_file,'file') || overwrite_coeff)
  %  for ll = 1:size(dyadic_spectrograms_reshaped,2)
  [COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(appearance_features_agg(:,:)');
  appearance_coeffs = COEFF;
   % end
    save(appearance_coeff_file,'appearance_coeffs');
else
load(appearance_coeff_file)
end
COEFF = appearance_coeffs;

centered_appearance = centered_appearance'*COEFF;
ML_features.appearance_features_agg_score_whitened = centered_appearance(:,1:num_eigenpcs);
clear appearance_features_agg centered_appearance

%% get selected differences
for ll = 1:numel(appearance_pairs)
   appearance_features(:,ll) = vectornorm( mocapstruct.markers_preproc.(appearance_pairs{ll}{1}),mocapstruct.markers_preproc.(appearance_pairs{ll}{2}),2);
end

ML_features.appearance_features = appearance_features(framelist_true,:);
clear appearance_features





%% trajectory features
% get average trunk displacement, velocity, acceleration, over different
% time windows: 10,30,100,300 frame windows
% also head head vel, arm and leg velocity
  clipped_pre_names = fieldnames(clipped_pre);

difforders = [10,33,100,300];
params.difforder = 10;
params.medfiltorder = 3;
params.gaussorder = 2.5;
ML_features.trunk_vel =zeros(numel(difforders),numel(framelist_true));
ML_features.head_vel =zeros(numel(difforders),numel(framelist_true));

velcomp_names = {'abs','x','y','z'};
absolute_velocity_names = {'trunk'};
absolute_velocity_markers = {[4:8]};

rel_velocity_names = {'head','trunk','hipL','hipR'};
rel_velocity_markers = {[1:3],[4,6:8],[9],[10]};
num_spectrogram_pcs= 15;
for ll = 1:numel(difforders)
    fprintf('starting absolute and relative velocity for windowsize %f \n',difforders(ll));
    
    params.difforder_movav = difforders(ll);
    for kk = 1:numel(absolute_velocity_names)
    [temp1,velcomp1,velstd,velstd_comp] = get_markersubgroup_velocity(mocapstruct.markers_preproc,absolute_velocity_markers{kk},params);
    ML_features.(strcat('absolute_velocity_',absolute_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = temp1(framelist_true);
    ML_features.(strcat('absolute_std_velocity_',absolute_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = velstd(framelist_true);

    for jj =1:size(velcomp1,2)
           ML_features.(strcat('absolute_velocity_',absolute_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =velcomp1(framelist_true,jj); 
             ML_features.(strcat('absolute_std_velocity',absolute_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =velstd_comp(framelist_true,jj); 

    end
    end
    
    
    %% relative velocity of the subgroups
    for kk = 1:numel(rel_velocity_names)
    [temp1,velcomp1,velstd,velstd_comp] = get_markersubgroup_velocity( clipped_pre,rel_velocity_markers{kk},params);
    ML_features.(strcat('rel_velocity_',rel_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = temp1(:);
        ML_features.(strcat('rel_std_velocity_',rel_velocity_names{kk},'_',velcomp_names{1},'_',num2str(difforders(ll)))) = velstd(:);

    for jj =1:size(velcomp1,2)
           ML_features.(strcat('rel_velocity_',rel_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =velcomp1(:,jj); 
               ML_features.(strcat('rel_std_velocity_',rel_velocity_names{kk},'_',velcomp_names{1+jj},'_',num2str(difforders(ll)))) =velstd_comp(:,jj); 

    end
    end
end

    %% spectrograms of the subgroups
  
%% load old coefficients or not
if (~exist(dyn_coeff_file_markers,'file') || overwrite_coeff)
        COEFFS_feat = cell(1,numel(rel_velocity_names));
else
load(dyn_coeff_file_markers)
end
    
%% get the pcs of the spectrogram of markers
    opts.clustering_overlap = 0;
    for kk = 1:numel(rel_velocity_names)
            fprintf('computing spectrograms for markers %f \n',kk)

       % should refactor this get_dyadic_pcs(clipped_pre,rel_velocity_markers{kk},opts)      
        agg_features_here = [];
        for ll = ( rel_velocity_markers{kk})
             agg_features_here = cat(2,agg_features_here,clipped_pre.(clipped_pre_names{ll}));
        end
             %  [dyadic_spectrograms,fr,~] = get_dyadic_spectrogram( agg_features_here',opts);
            dyadic_spectrograms = [];
             for jj = 1:size(agg_features_here,2)
                  [~,fr_temp,~,dyadic_spectrograms_temp] = spectrogram(agg_features_here(:,jj),opts.clustering_window,...
            opts.clustering_overlap,1:30,opts.fps);
        dyadic_spectrograms = cat(1,dyadic_spectrograms,dyadic_spectrograms_temp);
             end
              dyadic_spectrograms = log( dyadic_spectrograms);
                
              %% if need new coefficients -- these are constant across files
              if (~exist(dyn_coeff_file_markers,'file') || overwrite_coeff)
     [COEFFS_feat{kk}, dyadic_spectrograms_score, ~, ~,explained] = pca(squeeze(dyadic_spectrograms)');
              else
        dyadic_spectrograms_score = bsxfun(@minus,squeeze(dyadic_spectrograms),mean(squeeze(dyadic_spectrograms),2) )'*squeeze( COEFFS_feat{kk});
              end
              
replication_factor = floor(size(agg_features_here,1)./size(dyadic_spectrograms_score,1));
dynamics_pcs = repelem(dyadic_spectrograms_score(:,1:num_spectrogram_pcs),replication_factor,1);
dynamics_pcs = cat(1,dynamics_pcs,zeros(size(agg_features_here,1)-size(dynamics_pcs,1),size(dynamics_pcs,2)));

    ML_features.(strcat('spectrogram_pcs_',rel_velocity_names{kk})) = dynamics_pcs;
        ML_features.(strcat('spectrogram_coeffs_',rel_velocity_names{kk})) = COEFFS_feat{kk}(:,1:num_spectrogram_pcs);

    end
    
      if (~exist(dyn_coeff_file_markers,'file') || overwrite_coeff)
          save(dyn_coeff_file_markers,'COEFFS_feat')
      end
%     
%  figure(102)
%  imagesc(ML_features.spectrogram_coeffs_head)
%  
    
    
    %% wavelets of the subgroups
    do_wavelets = 0;
    if (do_wavelets)
        downsample = 3;
    maxframes = size(clipped_pre.HeadF,1);
    frames_use = 1:downsample: maxframes;
    clustering_ind = frames_use; %intersect with the chunking
    cluster_fps = 300./downsample;
        opts.num = 1; % number modes (spectrograms to find) (usually want full dimension)

%% setup cluster properties
opts.whiten = 0;
opts.frameNormalize = 0;
opts.clustermethod = 'GMM';
% num = pcuse;
opts.ds = 1; % down sampling
opts.samprate = 100;
opts.params = struct;
opts.params.samplingFreq = 100;
opts.params.numPeriods=25; %distinct number of frequencies to use
opts.params.minF = 0.3; % min freq to analyze
opts.params.maxF = 15; % max freq to analyze


% for gmm model parameters
opts.pcuse = 20;
opts.numclusters = 100;
opts.lambda = 0.1; % regularization

    
    tic
    [wavelets]= return_wavelets(clipped_pre.SpineL(clustering_ind,:),1:numel(clustering_ind),opts);
    fprintf('Time for wavelet cluster %f \n',toc)
    toc
    end
%     
% 
%     %relative velocities to mean
%         [temp2,velcomp2] = get_markersubgroup_velocity(mocapstruct.markers_aligned_preproc,[1:3],params);
%                 [temp3,velcomp3] = get_markersubgroup_velocity(mocapstruct.markers_aligned_preproc,[4:8],params);
%                 [temp4,velcomp4] = get_markersubgroup_velocity(mocapstruct.markers_aligned_preproc,[9],params);
%                 [temp5,velcomp5] = get_markersubgroup_velocity(mocapstruct.markers_aligned_preproc,[10],params);
%                % [temp6,velcomp6] = get_markersubgroup_velocity(mocapstruct.markers_aligned_preproc,[17],params);
% 
%       %  temp2 = get_markersubgroup_velocity(mocapstruct.markers_aligned_preproc,[1:3],params);
%        % temp2 = get_markersubgroup_velocity(mocapstruct.markers_aligned_preproc,[1:3],params);
% 
% ML_features.head_rel_vel(ll,:) = temp2(framelist_true);
% ML_features.trunk_rel_vel(ll,:) = temp3(framelist_true);
% ML_features.hipL_rel_vel(ll,:) = temp4(framelist_true);
% ML_features.hipR_rel_vel(ll,:) = temp5(framelist_true);

% 
% 
% ML_features.spectrogram_pcs_head(:,1:3)
%   
%    figure(98)
%    plot3(ML_features.spectrogram_pcs_trunk(1:100:end,10),ML_features.spectrogram_pcs_trunk(1:100:end,2),ML_features.spectrogram_pcs_trunk(1:100:end,3),'+')
% % plot3(dyadic_spectrograms_score(:,1),dyadic_spectrograms_score(:,2),dyadic_spectrograms_score(:,3),'+')
% % plot3(ML_features.spectrogram_pcs_hipR(:,1),ML_features.spectrogram_pcs_hipR(:,2),ML_features.spectrogram_pcs_hipR(:,3),'+')
% % 
% % 
%  subset = 1:100:100*7000;
%  mapped = tsne(cat(2,ML_features.spectrogram_pcs_head(subset,1:15),ML_features.spectrogram_pcs_trunk(subset,1:15),...
%    ML_features.spectrogram_pcs_hipL(subset,1:15),ML_features.spectrogram_pcs_hipR(subset,1:15)));
% 
% mapped = tsne(cat(2,ML_features.spectrogram_pcs_trunk(subset,1:15),...
%    ML_features.spectrogram_pcs_hipL(subset,1:15),ML_features.spectrogram_pcs_hipR(subset,1:15)));
% 
% figure(22)
% plot(mapped(:,1),mapped(:,2),'+')
% 
%      [COEFF_dyn, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(squeeze(dyadic_spectrograms)');
%      
%    
%       [dyadic_spectrograms,fr,timespect] = get_dyadic_spectrogram( clipped_pre.SpineL',opts);

% 
%  params.difforder_movav = 100;
% [vectorout,vectorstdout]= get_vector_velocity(clipped_pre2.HeadF(:,1),params);
% [vectorout2,vectorstdout]= get_vector_velocity(clipped_pre2.HeadF(:,3),params);
% 
% 
% figure(59)
% plot(10*vectorout,'b')
% hold on
% plot(10*vectorout2,'g')
% plot(10*ML_features.rel_velocity_head_abs_300,'k')
% 
% plot(clipped_pre.HeadF(:,3),'c')
% plot(clipped_pre.HeadB(:,3),'r')
% plot(clipped_pre.HeadL(:,3),'y')
% hold off
% 
% 
% figure(55)
% plot(ML_features.rel_velocity_head_abs_300,'g')
% hold on
% plot(ML_features.rel_velocity_hipL_abs_300,'b')
% plot(ML_features.rel_velocity_hipR_abs_300,'r')
% plot(ML_features.rel_velocity_trunk_abs_300,'k')
% hold off
% 
% plot(ML_features.rel_velocity_head_abs_33,'b')
% 
% plot(clipped_pre.HeadF(:,3),'c')
% plot(clipped_pre.HeadB(:,3),'r')
% plot(clipped_pre.HeadL(:,3),'y')
% hold off



fn_ml = fieldnames(ML_features);
abs_fn = intersect(intersect(find(cellfun(@numel,strfind(fn_ml,'100'))),find(cellfun(@numel,strfind(fn_ml,'rel_')))),find(cellfun(@numel,strfind(fn_ml,'abs'))));
abs_fn2 = intersect(find(cellfun(@numel,strfind(fn_ml,'100'))),find(cellfun(@numel,strfind(fn_ml,'_z_'))));
abs_fn3 = intersect(find(cellfun(@numel,strfind(fn_ml,'100'))),find(cellfun(@numel,strfind(fn_ml,'absolute_velocity_trunk_abs'))));

abs_fn = cat(1,abs_fn,abs_fn2,abs_fn3);
fn_ml(abs_fn)
agg_dyn = [];
for kk = 1:numel(abs_fn)
agg_dyn = cat(2,agg_dyn,ML_features.(fn_ml{abs_fn(kk)}));
end
agg_dyn = agg_dyn';
dyn_mean = mean(agg_dyn(:,:),2);
%center the pose
centered_dyn = bsxfun(@minus,agg_dyn,dyn_mean);
    agg_dyn(abs(agg_dyn)>20) = 20;
% if doesnt exist compute eigenposture feature coefficients and dynamics
% % coefficients
% if (~exist(dyn_coeff_file_2,'file') || overwrite_coeff)
        [COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(agg_dyn(:,:)');
%     save(pose_coeff_file,'COEFF');
% else
%     load(pose_coeff_file);
% end
%compute score
dyn_score = centered_dyn'*COEFF;
dyn_score_whitened = dyn_score;%bsxfun(@rdivide,pose_score,std(pose_score,[],1));
ML_features.dyn_score = dyn_score_whitened(:,1:5);%num_eigenpcs);

figure(44)
plot(ML_features.dyn_score(:,3))
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
        params.difforder_movav = difforders(ll);
   %   ML_features.pose_window_sd(ll,:,mm) = movstd(feature_mat_dyn_posespace_whitened(:,mm),difforders(ll));
        [ML_features.pose_velocity(ll,:,mm),ML_features.pose_window_sd(ll,:,mm)] = get_vector_velocity(feature_mat_dyn_posespace_whitened(:,mm),params);
       
    end
end

%% get joint angle features
%saggital/cross section
saggital_names = {'head_sagg','neck_sagg','spine_sagg'};
saggital_pairs =  {[2,3],[3,4],[4,5]}; %head, neck, spine angles , look in the z-y plane
%transverse/overhead
transverse_names = {'head_trans','neck_trans','spine_trans','hipl_trans','hipr_trans'};
transverse_pairs =  {[2,3],[3,4],[4,5],[5,6],[5,7]}; %head, neck, spine angles , look in the z-y plane
%coronal/along spine
coronal_names = {'head_coronal','hipl_coronal','hipr_coronal'};
coronal_pairs =  {[1,3],[5,6],[5,7]}; %head, neck, spine angles , look in the z-y plane


saggital_inds = [2,3];
coronal_inds = [1,3];
transverse_inds = [1,2];
%overall_inds = [1,2,3]; use on knees and arms
%transverse_pairs 

segment_pairs = {{'HeadB','HeadL'},{'HeadF','HeadB'},{'HeadB','SpineF'},{'SpineF','SpineM'} ,...
    {'SpineL','SpineM'},{'SpineL','HipL'},{'SpineL','HipR'},};    
    
    jointangle_struct = struct();
for ll = 1:numel(saggital_pairs)  
   vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(1)}{1})-...
       mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(1)}{2});
   
     vec2 =  mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(2)}{1})-...
       mocapstruct.markers_aligned_preproc.(segment_pairs{saggital_pairs{ll}(2)}{2});
  % jointangle_struct.saggital_names{ll} = [];
   jointangle_struct.(saggital_names{ll}) = ...    
       acos(dot(vec1(framelist_true,saggital_inds)',vec2(framelist_true,saggital_inds)')...
       ./(sqrt(sum(vec1(framelist_true,saggital_inds).^2,2)).*sqrt(sum(vec2(framelist_true,saggital_inds).^2,2)))');
end

for ll = 1:numel(coronal_pairs)  
   vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(1)}{1})-...
       mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(1)}{2});
   
     vec2 =  mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(2)}{1})-...
       mocapstruct.markers_aligned_preproc.(segment_pairs{coronal_pairs{ll}(2)}{2});
  % jointangle_struct.saggital_names{ll} = [];
   jointangle_struct.(coronal_names{ll}) = ...    
       acosd(dot(vec1(framelist_true,saggital_inds)',vec2(framelist_true,coronal_inds)')...
       ./(sqrt(sum(vec1(framelist_true,coronal_inds).^2,2)).*sqrt(sum(vec2(framelist_true,coronal_inds).^2,2)))');
end


for ll = 1:numel(transverse_pairs)  
   vec1 =  mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{1})-...
       mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(1)}{2});
   
     vec2 =  mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(2)}{1})-...
       mocapstruct.markers_aligned_preproc.(segment_pairs{transverse_pairs{ll}(2)}{2});
  % jointangle_struct.saggital_names{ll} = [];
   jointangle_struct.(transverse_names{ll}) = ...    
       acosd(dot(vec1(framelist_true,transverse_inds)',vec2(framelist_true,transverse_inds)')...
       ./(sqrt(sum(vec1(framelist_true,transverse_inds).^2,2)).*sqrt(sum(vec2(framelist_true,transverse_inds).^2,2)))');
end

ML_features.joint_angles = jointangle_struct;

figure(66)
plot(( jointangle_struct.(transverse_names{1})))




%% pcs here


%% dynamics of the joint angles -- are they cleaner?
head_angles = {'head_sagg','head_trans','neck_sagg','neck_trans','head_coronal'};
trunk_angles = {'spine_sagg','hipl_coronal','hipr_coronal','spine_trans','hipl_trans','hipr_trans'};

angle_lists = {head_angles,trunk_angles};
angle_list_name = {'head_angle','trunk_angle'};
ML_features.trunk_angles = trunk_angles;


%% load old coefficients or not
if (~exist(jadyn_coeff_file,'file') || overwrite_coeff)
        COEFFS_feat = cell(1,numel(angle_lists));
else
load(jadyn_coeff_file)
end

  opts.clustering_overlap = 75;
 for kk = 1:numel(angle_lists)
     ML_features.angle_names{kk} = (angle_lists{kk});
            fprintf('computing spectrograms for markers %f \n',kk)

       % should refactor this get_dyadic_pcs(clipped_pre,rel_velocity_markers{kk},opts)      
        agg_features_here = [];
        for ll = 1:numel( angle_lists{kk})
             agg_features_here = cat(2,agg_features_here,ML_features.joint_angles.( angle_lists{kk}{ll})');
        end
             %  [dyadic_spectrograms,fr,~] = get_dyadic_spectrogram( agg_features_here',opts);
            dyadic_spectrograms = [];
             for jj = 1:size(agg_features_here,2)
          %        [~,fr_temp,~,dyadic_spectrograms_temp] = spectrogram(agg_features_here(:,jj),opts.clustering_window,...
         %   opts.clustering_overlap,1:30,opts.fps);
        
        params.tapers = [5 7];
        params.Fs = 300;
        
        [S,t,f]=mtspecgramc( bsxfun(@minus,agg_features_here(:,jj),mean(agg_features_here(:,jj))),[0.5 0.5],params);
        dyadic_spectrograms_temp = (S(:,1:50).^2)';
        ML_features.mean_ja_spect{kk}(jj,:) = log(mean(dyadic_spectrograms_temp,2)');
        
%      figure(33)
%         imagesc(t,f,log10(S.^2)')
%         ylim([0 30])
%         %S2 = dyadic_spectrograms_temp';
%         figure(34)
%         induse = 1:30;
%         plot(f(induse),mean(log10(S(:,1:30)),1))
%         
        dyadic_spectrograms = cat(1,dyadic_spectrograms,dyadic_spectrograms_temp);
             end
              dyadic_spectrograms = log( dyadic_spectrograms);
         
                %% if need new coefficients -- these are constant across files
              if (~exist(jadyn_coeff_file,'file') || overwrite_coeff)
     [COEFFS_feat{kk}, dyadic_spectrograms_score, ~, ~,explained] = pca(squeeze(dyadic_spectrograms)');
             %dyadic_spectrograms_score = bsxfun(@minus,squeeze(dyadic_spectrograms)',mean(squeeze(dyadic_spectrograms)',2) )*squeeze( COEFFS_feat{kk});

              else
        dyadic_spectrograms_score = bsxfun(@minus,squeeze(dyadic_spectrograms),mean(squeeze(dyadic_spectrograms),2) )'*squeeze( COEFFS_feat{kk});
              end
              
    % [COEFFS_feat, dyadic_spectrograms_score, ~, ~,explained] = pca(squeeze(dyadic_spectrograms)');
     
replication_factor = floor(size(agg_features_here,1)./size(dyadic_spectrograms_score,1));
dynamics_pcs = repelem(dyadic_spectrograms_score(:,1:num_spectrogram_pcs),replication_factor,1);
dynamics_pcs = cat(1,dynamics_pcs,zeros(size(agg_features_here,1)-size(dynamics_pcs,1),size(dynamics_pcs,2)));

 ML_features.ja_freq = f(1:50);
   ML_features.(strcat('spectrogram_pcs_',angle_list_name{kk})) = dynamics_pcs;
      ML_features.(strcat('spectrogram_coeffs_',angle_list_name{kk})) = COEFFS_feat{kk}(:,1:num_spectrogram_pcs);

 end

 
num_pc = 6;
num_angle = 5;

 for pc_plot = 1:6
 kk_plot = 1;
 
COEFFS_resh= reshape(COEFFS_feat{kk_plot}(:,pc_plot),numel( ML_features.ja_freq),[]);
deviation = std(ML_features.(strcat('spectrogram_pcs_',angle_list_name{kk_plot}))(:,pc_plot),[],1);

 summed_coeffs =  ML_features.mean_ja_spect{kk_plot}'+COEFFS_resh*deviation;
  summed_coeffs_minus =  ML_features.mean_ja_spect{kk_plot}'-COEFFS_resh*deviation;

   summed_coeffs_exp =  (10.^(ML_features.mean_ja_spect{kk_plot}'+COEFFS_resh*deviation)-10.^(ML_features.mean_ja_spect{kk_plot}'));

  
  for angle_plot = 1:5
  
  
 figure(44)
 subplot(num_angle,num_pc,pc_plot+num_pc*(angle_plot-1))
 plot(COEFFS_resh(:,angle_plot)*deviation,'r');
 hold on
 %plot(ML_features.mean_ja_spect{kk_plot}(angle_plot,:)','k');
 plot(-COEFFS_resh(:,angle_plot)*deviation,'b' )

 if (pc_plot == 1)
     ylabel(ML_features.angle_names{kk_plot}{angle_plot})
 end
 if (angle_plot == 1)
     ntitle(strcat('PC ',num2str(pc_plot)));
 end
 
 
 
 figure(45)
 subplot(num_angle,num_pc,pc_plot+num_pc*(angle_plot-1))
 plot(summed_coeffs(:,angle_plot),'r');
 hold on
 plot(ML_features.mean_ja_spect{kk_plot}(angle_plot,:)','k');
 plot(summed_coeffs_minus(:,angle_plot),'b' )
 
 if (pc_plot == 1)
     ylabel(ML_features.angle_names{kk_plot}{angle_plot})
 end
 if (angle_plot == 1)
     ntitle(strcat('PC ',num2str(pc_plot)));
 end
 
 figure(46)
 subplot(num_angle,num_pc,pc_plot+num_pc*(angle_plot-1))
 plot(exp(summed_coeffs(:,angle_plot)),'r');
 hold on
  plot(exp(ML_features.mean_ja_spect{kk_plot}(angle_plot,:)'),'k');
  plot(exp(summed_coeffs_minus(:,angle_plot)),'b');

 if (pc_plot == 1)
     ylabel(ML_features.angle_names{kk_plot}{angle_plot})
 end
 if (angle_plot == 1)
     ntitle(strcat('PC ',num2str(pc_plot)));
 end
  end
 end
      if (~exist( jadyn_coeff_file,'file') || overwrite_coeff)
          save( jadyn_coeff_file,'COEFFS_feat')
      end


    
    figure(98)
   plot3(ML_features.spectrogram_pcs_head_angle(1:100:100*10000,1),ML_features.spectrogram_pcs_head_angle(1:100:100*10000,2),ML_features.spectrogram_pcs_head_angle(1:100:100*10000,6),'+')
% % % plot3(dyadic_spectrograms_score(:,1),dyadic_spectrograms_score(:,2),dyadic_spectrograms_score(:,3),'+')
% % % plot3(ML_features.spectrogram_pcs_hipR(:,1),ML_features.spectrogram_pcs_hipR(:,2),ML_features.spectrogram_pcs_hipR(:,3),'+')
% % % 
% % % 
%   subset = 1:100:100*2400;
%  % mapped = tsne(cat(2,ML_features.spectrogram_pcs_head(subset,1:15)));   
%    mapped = tsne(cat(2,ML_features.spectrogram_pcs_head(subset,1:15),ML_features.spectrogram_pcs_trunk(subset,1:15)));   
% 
%    figure(22)
% plot(mapped(:,1),mapped(:,2),'+')
%   
%  figure(102)
%  imagesc(ML_features.spectrogram_coeffs_head)
%  


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