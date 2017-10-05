% OLD OPTIONS FOR VISUALIZING DIFFERENT DATASETS
% save_tag = 'Vicon8pre_task_prelesion';
% plotdirectory = strcat(savedirectory,save_tag,filesep);
% mkdir(plotdirectory)
% 
% %create second as well
% save_tag = 'Vicon8_caff';
% plotdirectory2 = strcat(savedirectory,save_tag,filesep);
% mkdir(plotdirectory2)
% 
% %create second as well
% save_tag = 'Vicon8_amph';
% plotdirectory3 = strcat(savedirectory,save_tag,filesep);
% mkdir(plotdirectory3)
% 
% save_tag = 'Vicon8_dlslesion_early';
% plotdirectory4 = strcat(savedirectory,save_tag,filesep);
% mkdir(plotdirectory4)
% 
% save_tag = 'Vicon8_veh';
% plotdirectory5 = strcat(savedirectory,save_tag,filesep);
% mkdir(plotdirectory5)
% 
% save_tag = 'JDM25pre';
% plotdirectory6 = strcat(savedirectory,save_tag,filesep);
% mkdir(plotdirectory6)
% 
% 
% 
% goodinds_2 = find(cellfun(@numel,strfind(mocapfilestruct.PreLesion.mocapfiles{6},'caff')));
% mocapfilearray_2 = mocapfilestruct.PreLesion.mocapfiles{6}(goodinds_2);
% 
% goodinds_3 = find(cellfun(@numel,strfind(mocapfilestruct.PreLesion.mocapfiles{7},'amph')));
% mocapfilearray_3 = mocapfilestruct.PreLesion.mocapfiles{7}(goodinds_3);
% 
% goodinds_4 = find(cellfun(@numel,strfind(mocapfilestruct.UniLesion.mocapfiles{2},'overnight_two1')));
% mocapfilearray_4 = mocapfilestruct.UniLesion.mocapfiles{2}(goodinds_4);
% 
% goodinds_5 = find(cellfun(@numel,strfind(mocapfilestruct.PreLesion.mocapfiles{5},'veh')));
% mocapfilearray_5 = mocapfilestruct.PreLesion.mocapfiles{5}(goodinds_5);
% 
% goodinds_6 = find(cellfun(@numel,strfind(mocapfilestruct.PreLesion.mocapfiles{6},'overnight1')));
% mocapfilearray_6 = mocapfilestruct.PreLesion.mocapfiles{6}(goodinds_6);
% % 
% % goodinds_2 = find(cellfun(@numel,strfind(mocapfilestruct.UniLesion.mocapfiles{5},'overnight1')));
% % mocapfilearray_2 = mocapfilestruct.UniLesion.mocapfiles{5}(goodinds_2);
% 
% [mocapstruct_2] = preprocess_mocap_data(mocapfilearray_2,mocapfilestruct.mocapdir);
% [mocapstruct_3] = preprocess_mocap_data(mocapfilearray_3,mocapfilestruct.mocapdir);
% [mocapstruct_4] = preprocess_mocap_data(mocapfilearray_4,mocapfilestruct.mocapdir);
% [mocapstruct_6] = preprocess_mocap_data(mocapfilearray_6,mocapfilestruct.mocapdir);


%% plotting directory -- change this
%mocapmasterdirectory = 'E:\Bence\Data\Motionanalysis_captures\';
mocapmasterdirectory = '\\140.247.178.37\Jesse\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots',filesep);
mkdir(savedirectory);

%% load or create struct
%createmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);

%% get the desired files
descriptor_struct_1 = struct();
descriptor_struct_1.day = 5;
descriptor_struct_1.tag = 'overnight1';
descriptor_struct_1.cond = 'PreLesion';

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);

%% either load  or preprocess from scratch
[mocapstruct] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1);

%% can also load from saved data
%[mocapstruct_load] = load_mocap_data(mocapfilearray,mocapfilestruct.mocapdir);
comp_fraction_moving(mocapstruct)

% plot fraction time missing
plotfractionmissing(mocapstruct)

plot_marker_characteristics(mocapstruct)
%% code doesn't exist yet
%do_task_analysis(mocapstruct)
%compare_task(mocapstruct,mocapstruct_2);



%% EXAMPLE to visualize a simple portion of the trace
% traceuse = mocapstruct.markers_preproc.HeadF(:,3);
   params_clip.fps = mocapstruct.fps;
   [markers_clipped,clipped_index_here] = hipass_clip(mocapstruct.markers_preproc,mocapstruct.bad_frames_agg{1},params_clip);
traceuse = markers_clipped.HeadF(:,2);
%regionuse = 170000:174000;
figure(111)
plot(0:1./300:(numel(traceuse)-1)./300,traceuse)
xlabel('Time (s)')
ylabel('Head Z (mm)')
%print('-dpng',strcat(plotdirectory,filesep,'head_vel_ex.png'))

% do_simple_analysis(traceuse(regionuse),veluse(regionuse),clipped_index_here(regionuse),'head',fps,markers_preproc,plotdirectory,markercolor,links)

%% compute 


%% do clustering
cluster_here = [2];
    downsample = 3;
savedirectory_subcluster = strcat(mocapstruct.plotdirectory,filesep,'subclusterplots_',num2str(cluster_here),filesep);
mkdir(savedirectory_subcluster);

[modular_cluster_properties] = get_modularclustproperties(mocapstruct);

%% the clipped index is the index of clipped frames in the full trace
[modular_cluster_properties] = get_clustering_features(mocapstruct,modular_cluster_properties,2)  ;

    %% first compute the eigenposes to look for postural differences
      %% look at the eigenposes of the data
        fprintf('computing eigenposes \n')
  
  eigenpose = compute_eigenposes(mocapstruct,modular_cluster_properties,cluster_here);
    
   %% GMM cluster
   opts.clustering_window = mocapstruct.fps;
   opts.clustering_overlap = floor(mocapstruct.fps*0.5);
   opts.fps = mocapstruct.fps;
   opts.num_pcs_1 = 50;
   opts.num_pcs_2 = 30;
   opts.num_clusters = 200;
   
   %% get clusters and metrics
[cluster_struct_spect] = Cluster_GMM(modular_cluster_properties.agg_features{cluster_here},opts,1:size(modular_cluster_properties.agg_features{cluster_here},2 ));
cluster_struct_spect.feature_labels = modular_cluster_properties.feature_labels{cluster_here};
cluster_struct_spect.clipped_index_agg = modular_cluster_properties.clipped_index{cluster_here};               
 % [cluster_struct_spect] = clusterMetricsTodd_JDM(cluster_struct_spect);
  
  %have to convert temporally 
  %cluster_struct_GMM=convert_cluster_GMM_struct(cluster_struct_GMM)
  
  do_movies=1;
do_cluster_plots=1;
               plot_cluster_means_movies(savedirectory_subcluster,cluster_struct_spect,modular_cluster_properties,cluster_here,...
mocapstruct.markers_preproc,do_movies,mocapstruct,do_cluster_plots) 
  %% have to add other spectrogram features to get plots/movies
  make_cluster_descriptors(cluster_struct_spect,modular_cluster_properties.agg_features{cluster_here}(3,:),savedirectory_subcluster)

  
figure(222)
[n,x] = hist((cluster_struct_spect.labels),opts.num_clusters);
bar(x,sort(n)./sum(n))
set(gca,'YScale','log')
ylabel('Cluster Frequency')
xlabel('Cluster Number')
box off
  print('-dpng',strcat(mocapstruct.plotdirectory,'cluster_frequency_spectrogram.png')) 
    


  %% wavelet cluster and clusterobj

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
opts.params.maxF = 20; % max freq to analyze

% for gmm model parameters
opts.pcuse = 20;
opts.num = 100;
opts.numclusters = opts.num;
opts.lambda = 0.1; % regularization
clustering_ind = 1:downsample:size(modular_cluster_properties.agg_features{2},2);

cluster_struct_wavelet = WaveletCluster(modular_cluster_properties.agg_features{2}(:,clustering_ind)',1:numel(clustering_ind),opts);    
              cluster_struct_wavelet.feature_labels = modular_cluster_properties.feature_labels{cluster_here};
                cluster_struct_wavelet.clipped_index = modular_cluster_properties.clipped_index{cluster_here};
      [cluster_struct_wavelet] = clusterMetricsTodd_JDM(cluster_struct_wavelet);
      
do_movies=1;
do_cluster_plots=0;
                plot_cluster_means_movies(savedirectory_subcluster,cluster_struct_wavelet,modular_cluster_properties.subcluster_names{cluster_here},...
                    mocapstruct.markers_preproc,do_movies,mocapstruct,do_cluster_plots) 

  make_cluster_descriptors(cluster_struct_wavelet,modular_cluster_properties.agg_features{cluster_here}(3,:),savedirectory_subcluster)
  cluster_scan = [30,50];
  % compare across a few cluster numbers
   cluster_struct_wavelet.hyperparam_comparison_struct = compare_numberofclusters(cluster_struct_wavelet,modular_cluster_properties,cluster_here,  cluster_scan,opts);
   %agg_features,cluster_scan)
%  cluster_struct = makeclusterstruct(struct(),labels,feature_mu,feature_sigma,fr,feat_pcs,number_of_clust(mmm),...
%         feature_labels,clustering_ind,clipped_index_agg{mmm});
    

            
  
  %% compare two clusters
       %  feature_euclidean_distance(cluster_struct,cluster_struct2,wtAll,wtAll2,savedirectory_subcluster)
       
       
       
       
       %% compute the cluster compactness (good in principle)
          % concentration_index= zeros(2,num_clusters);
    %
    % for zz = 1:num_clusters
    %    concentration_index(1,zz) = sum((squeeze(sum(abs(feature_mu_reshaped(zz,:,:)),2))./sum(squeeze(sum(abs(feature_mu_reshaped(zz,:,:)),2)))).^2)...
    %     ./sum(squeeze(sum(abs(feature_mu_reshaped(zz,:,:)),2)));
    %
    %
    %    concentration_index(2,zz) = sum(sum(abs(feature_mu_reshaped(zz,:,:))./sum(squeeze(sum(abs(feature_mu_reshaped(zz,:,:)),2))).^2))...
    %     ./sum(squeeze(sum(abs(feature_mu_reshaped(zz,:,:)),2)));
    %
    % end
    %
    % title('weighting concentration')
    % figure(344)
    % hist(squeeze(concentration_index(2,:)),20)
    %
    % for zz = 1:num_clusters
    %     cluster_frames = clustering_ind(find(labels == zz));
    %     %% first compute pose similarity
    % for jj = 1:numel(cluster_marker_inds)
    %    % average_post_struct(zz).(marker_names{jj}) = nanmean(markers_preproc.(marker_names{ll}))
    %    % have to compute from centered/aligned clusters
    %     for kk = 1:numel(cluster_marker_inds)
    % % average_pose(1,jj,kk,zz) = nanmean(squeeze(delta_markers(cluster_marker_inds(jj),cluster_marker_inds(kk),cluster_frames,4)));
    % % average_pose(2,jj,kk,zz) = nanstd(squeeze(delta_markers(cluster_marker_inds(jj),cluster_marker_inds(kk),cluster_frames,4)));
    % %average_pose(3,jj,kk,zz) = nanstd(squeeze(delta_markers(jj,kk,:,4)));
    %
    %     end
    % end
