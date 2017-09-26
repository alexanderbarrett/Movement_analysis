% read_c3d_file_better
set(0,'DefaultFigureWindowStyle','docked')

%% Vicon8
mocapdir = 'E:\Bence\Data\Motionanalysis_captures\Vicon8\';
mocapnames = {'20170816\Generated_C3D_files\','20170817\Generated_C3D_files\','20170818\Generated_C3D_files\',...
    '20170820\Generated_C3D_files\',...
    '20170821\Generated_C3D_files\','20170822\Generated_C3D_files\','20170823\Generated_C3D_files\',... day7
    '20170824\Generated_C3D_files\','20170825\Generated_C3D_files\',...%8,9
    '20170827\Generated_C3D_files\','20170828\Generated_C3D_files\',...%10,11
    '20170830\Generated_C3D_files\','20170831\Generated_C3D_files\','20170901\Generated_C3D_files\','20170904\Generated_C3D_files\'};
filepath_array = [];
filepath_ind = 1;

days_to_compare{1} = 13;%[10:14];
%days_to_compare{2} = [10];

lever_thresholded_agg = [];
resample_analog_agg = [];
analog_agg = struct();
marker_agg = struct();

%task pre lesion: day 5, '*nolj*Recording*overnight5*.c3d'
%behavior pre lesion '*nolj*Recording_day7_caffeine3*.c3d'

for jj = days_to_compare{1}%9:numel(mocapnames)
   % filenames_here = dir(strcat(mocapdir,mocapnames{jj},'*nolj*Recording_day7_caffeine3*.c3d'));
        filenames_here = dir(strcat(mocapdir,mocapnames{jj},'*nolj*recording*overnight1*.c3d'));

    
    numel(filenames_here)
    for kk = 1:numel(filenames_here)
        if (numel(strfind(filenames_here(kk).name,'nolj')))
            filepath_array{filepath_ind} = strcat(mocapdir,mocapnames{jj},filenames_here(kk).name);
            filepath_ind = filepath_ind+1;
        end
    end
end

savedirectory = 'E:\Bence\Data\MOCAP\regression_animations\';
save_tag = 'Vicon8_task_postlesion';
plotdirectory = strcat(savedirectory,save_tag,filesep);
mkdir(plotdirectory)

[mocap_datecreate,sorted_mocap_serialtimes,filepath_array_sorted ] = sort_mocap_files(filepath_array,mocapdir);

% btkCloseAcquisition(acq);
fps = 300;

analog_factor = 1;
analog_fps = fps*analog_factor;

desired_length = 8*10^6;
chunksize = 300*fps;

%filepath_array_sorted_analog = get_analogactive_files(filepath_array_sorted);

[markers,analog,resample_analog,lever_thresholded,filestartpts] = concatenate_andsample_c3d(filepath_array_sorted,fps,analog_fps,...
    desired_length,chunksize);


%% get basic information about the dataset
%plot_marker_lever_xcorr(markers,lever_thresholded,analog);


%possible options
[markercolor,links] = day_loading_header('vicon8_20marker');

%% touch up the arms so that the elbow to arm is pointed away from the body
%because of elbow and arm swaps, also set to zero if only one of the limbs
%is observed
[markers,markers_aligned] = align_hands_elbows(markers,fps);



%% Data cleaning
% get frames to analyze based on set conditions
marker_names = fieldnames(markers);
marker_frame_length = size(markers.(marker_names{1}),1);
markers_preproc = markers;
markers_preproc_aligned = markers_aligned;


%% print the 'up times' for markers, and get the missing frames
% get rid of frames where markers are absent
missing_times = cell(1,numel(marker_names));
missing_times_postpreprocess = cell(1,numel(marker_names));

for ll = 1:numel(marker_names)
    missing_times{ll} = find(markers.(marker_names{ll})(:,1)==0);
    fprintf('For marker %s percentage of frames present %e \n',marker_names{ll},...
        100.*numel(missing_times{ll})./marker_frame_length);
end

fprintf('Percentage where two or more are gone %e \n',0)
bad_times = cat(1,missing_times{:});
good_times = setxor(1:marker_frame_length,bad_times);


fakedata = struct();
for fn = fieldnames(markers_preproc)'
    fakedata.(fn{1}) = cell(1,1);
end


%% other interpolation etc.
for ll = 1:numel(marker_names)
    fprintf('starting interpolation for markers, over 30 frames (100 ms) maximum %f \n',ll)
    % for jkj = 1:3 %need to do over all x y z simul and add spike correction
    [temp,fake_frames] = markershortinterp((markers_preproc.(marker_names{ll})),30,5);
    fakedata.(marker_names{ll}){1} = fake_frames;
    % if numel(temp)
    markers_preproc.(marker_names{ll}) = temp;
    %  end
    clear fake_frames
    
    fprintf('interpolating aligned markers \n')
    [temp,~] =markershortinterp((markers_preproc_aligned.(marker_names{ll})),100,5);
    % if numel(temp)
    markers_preproc_aligned.(marker_names{ll}) = temp;
    %  end
end

%% median filter the data to remove spikes
fprintf('median filtering %f \n')

for ll = 1:numel(marker_names)
    markers_preproc.(marker_names{ll}) = medfilt2(markers_preproc.(marker_names{ll}),[3,1]);
end


for ll = 1:numel(marker_names)
    missing_times_postpreprocess{ll} = find(markers_preproc.(marker_names{ll})(:,1)==0);
    fprintf('For marker %s percentage of missing frames %e \n',...
        marker_names{ll},numel(missing_times_postpreprocess{ll})./marker_frame_length);
end

figure(555)
bar(cellfun(@numel,missing_times_postpreprocess)./marker_frame_length)
box off
set(gca,'XTick',1:numel(marker_names),'XTickLabels',marker_names)
rotateXLabels(gca,45)
xlabel('Marker')
ylabel('Fraction of time tracked')
ylim([0 0.5])
xlim([0 numel(marker_names)+1])
print('-dpng',strcat(plotdirectory,'MarkerTracking.png'))
%set(get(gca,'XTickLabel'),'Rotation',-45);



%% Transform the data into relevant features for clustering and visualization

%'big-data' features
% get relative marker positions to one another (x,y,z)
num_markers = numel(markers);
marker_velocity = zeros(num_markers,marker_frame_length,4);
marker_position = zeros(num_markers,marker_frame_length,3);
abs_velocity_antialiased = zeros(num_markers,marker_frame_length);


dH = designfilt('lowpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 60/(fps/2), ...
    'DesignMethod', 'butter');
[f1,f2] = tf(dH);

%delta_markers_reshaped = [];
fprintf('getting velocities \n')
for ll = 1:numel(marker_names)
    marker_position(ll,:,1:3) = markers_preproc.(marker_names{ll});
    for mk = 1:3
    end   
    marker_velocity(ll,2:(end),1:3) = diff(markers_preproc.(marker_names{ll}));
    marker_velocity(ll,1,1:3) = marker_velocity(ll,2,1:3);
    marker_velocity(ll,:,4) = sqrt(sum((squeeze( marker_velocity(ll,:,1:3))).^2,2));
    abs_velocity_antialiased(ll,:) =  filtfilt(f1,f2, marker_velocity(ll,:,4));

    for lk = (ll+1):num_markers
        distance_here =   (markers_preproc.(marker_names{ll})-markers_preproc.(marker_names{lk}));
    end
end


%get aggregate feature matrix
%% simple bad frame detector
fprintf('finding bad frames \n')
bad_frames_agg = getbadframes(marker_velocity,marker_position,fps);
clear marker_position

%% rest move discrim
hist(conv(abs_velocity_antialiased(5,:),ones(1,300)./300),500)
% plot(conv(abs_velocity_antialiased(5,:),ones(1,300)./300))
% vel_thresh = 0.01;


%% visualize the difference between the ro
%bad_frames_elbowspine = unique(cat(2,bad_frames_agg{4},bad_frames_agg{5},bad_frames_agg{11},bad_frames_agg{12}));
%goodframeshere =setxor(1:size(markers_preproc.ArmL,1),bad_frames_elbowspine);


%mocapfilestruct = loadmocapfilestruct;
%preprocess data
%[bad_frames,move_frames,markers,markers_preproc,markers_aligned_preproc] = load_preprocessed_data(day,condition,nframes);


%% do a simple feature plot
traceuse = markers_preproc.HeadF(:,3);
    params_clip.fps = fps;
    [markers_clipped,clipped_index_here] = hipass_clip(markers_preproc,bad_frames_agg{1},params_clip);
veluse = markers_clipped.HeadF(:,3);
regionuse = 170000:174000;
traceuse = traceuse(clipped_index_here);
figure(111)
plot(traceuse)
do_simple_analysis(traceuse(regionuse),veluse(regionuse),clipped_index_here(regionuse),'head',fps,markers_preproc,plotdirectory,markercolor,links)


%% do task analysis, possibly save to plotdirectory
fprintf('doing task analysis \n')
params_task.fps = fps;
do_task_analysis(analog,lever_thresholded,params_task,markers_preproc,marker_velocity,plotdirectory,marker_names,markercolor,links)



%% plot marker characteristics to get a feel for overall velocities, range of motion, etc. 
fprintf('plotting marker characteristics \n')
params.fps= fps;
%find_rest_frames = (markers_preproc,marker_velocity,
plot_marker_characteristics(markers_preproc,marker_velocity,bad_frames_agg,params);


cluster_here = [2,7,1];
cluster_objects = cell(1,number_of_subclusters);
[cluster_markersets,subcluster_names,badframe_inds_agg,clustering_inds_agg] = get_modularclustproperties(mocapstruct)

fprintf('starting clustering \n')

%% Reduce dimensionality ,compute feature spectrogram, cluster
do_clustering_analysis = 1;


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
opts.numclusters = 100;
opts.lambda = 0.1; % regularization



%% store feature weights
wtAll_agg = cell(1,number_of_subclusters);
labels_agg = cell(1,number_of_subclusters);
feature_mu_agg = cell(1,number_of_subclusters);
clipped_index_agg = cell(1,number_of_subclusters);

feature_labels_agg = cell(1,number_of_subclusters);
clustering_ind_agg = cell(1,number_of_subclusters);
for mmm=2;%cluster_here;%[1:number_of_subclusters]%:numel(cluster_markersets)]
    
    fprintf('For subcluster %f Fraction of frames to exclude %f \n',mmm,(trace_length-numel(clustering_inds_agg{mmm}))./trace_length);
    
    % make plot folders
    savedirectory_subcluster =strcat(savedirectory, save_tag,num2str(mmm),filesep);
        savedirectory_subcluster2 =strcat(savedirectory, save_tag,num2str(mmm),'secondinterval',filesep);

    mkdir(savedirectory_subcluster)
    
    cluster_marker_inds = cluster_markersets{mmm};
    num_subcluster_markers = numel(cluster_marker_inds);
    marker_names_subcluster = marker_names(cluster_marker_inds);

    %% clip the makrers
    fprintf('clipping and hipass filtering marker positions \n')
    params_clip.fps = fps;
    [markers_clipped,clipped_index_agg{mmm}] = hipass_clip(markers_preproc_aligned,badframe_inds_agg{mmm},params_clip);
    
    agg_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));
    
    for ll = cluster_marker_inds      
        agg_features = cat(1,agg_features,markers_clipped.(marker_names{ll})');
    end
    
    
    
    %% change labels for features
    feature_labels = fieldnames(markers);
    feature_labels = feature_labels((cluster_marker_inds));
    feature_labels_agg{mmm} = feature_labels;
    
    %% downsample for clustering
    downsample = 3;
    maxframes = size(agg_features,2);
    frames_use = 1:downsample:size(agg_features,2);
    clustering_ind = intersect(0.1*300*10^4:0.8*300*10^4,frames_use); %intersect with the chunking
    clustering_ind_2 = intersect(0.8*300*10^4:1.4*300*10^4,frames_use); %intersect with the chunking
    
    cluster_fps = fps./downsample;
    clustering_ind_agg{mmm}= clustering_ind;
    
    opts.numclusters =   number_of_clust(mmm);
    opts.num = size(agg_features,1); % number modes (spectrograms to find) (usually want full dimension)
    
    
    %% look at the eigenposes of the data
        fprintf('computing eigenposes \n')

    agg_pose_features = [];%cat(1,squeeze(marker_velocity(:,1:frame_length_here,4)));
    
    for ll = cluster_marker_inds      
        agg_pose_features = cat(1,agg_pose_features,markers_preproc_aligned.(marker_names{ll})(clipped_index_agg{mmm},:)');       
    end
        
  eigenpose = computer_eigenposes(agg_pose_features,clustering_ind_2);
  
    
    
    
        %% compute the clusters
            fprintf('computing clusters \n')

            %% first do a hyperparameter test
            
    %% then cluster
    tic
    [labels , feature_mu, feature_sigma,wtAll,fr,score,feat_pcs]= WaveletCluster(agg_features(:,clustering_ind)',1:numel(clustering_ind),opts);
    fprintf('Time for wavelet cluster %f \n',toc)
    
    [labels2 , feature_mu2, feature_sigma2,wtAll2,fr2,score2,feat_pcs2]= WaveletCluster(agg_features(:,clustering_ind_2)',1:numel(clustering_ind_2),opts);
    fprintf('Time for wavelet cluster %f \n',toc)
    
    
    %% make cluster structs
    cluster_struct = makeclusterstruct(struct(),labels,feature_mu,feature_sigma,fr,feat_pcs,number_of_clust(mmm),...
        feature_labels,clustering_ind,clipped_index_agg{mmm});
        
     cluster_struct2 = makeclusterstruct(struct(),labels2,feature_mu2,feature_sigma2,fr2,...
         feat_pcs2,number_of_clust(mmm),feature_labels,clustering_ind_2,clipped_index_agg{mmm});
   
    
    %% make plots and get hyperparameters for cluster 1
      make_cluster_descriptors(cluster_struct,wtAll,clustering_ind,agg_features(3,:),savedirectory_subcluster,fps)
  
      [cluster_struct.shortCount_agg, cluster_struct.markovLLR_agg,...
                cluster_struct.entropy_agg, cluster_struct.exits_agg,...
                cluster_struct.meanDwell_agg, cluster_struct.meanClusterDist_agg] = clusterMetricsTodd_JDM(labels, wtAll);
                plot_cluster_means_movies(savedirectory_subcluster,cluster_struct,agg_features,wtAll,markers_preproc,do_movies) 

  make_cluster_descriptors(cluster_struct2,wtAll2,clustering_ind_2,agg_features(3,:),savedirectory_subcluster)
    
    
             [cluster_struct2.shortCount_agg, cluster_struct2.markovLLR_agg,...
                cluster_struct2.entropy_agg, cluster_struct2.exits_agg,...
                cluster_struct2.meanDwell_agg, cluster_struct2.meanClusterDist_agg] = clusterMetricsTodd_JDM(labels2, wtAll2);
            
              feature_euclidean_distance(cluster_struct,cluster_struct2,wtAll,wtAll2,savedirectory_subcluster)

                plot_cluster_means_movies(savedirectory_subcluster2,cluster_struct2,agg_features,wtAll2,markers_preproc,do_movies) 

                
            
            %% look at the feature euclidean distances
            
              %% store wavelet fetures            
    wtAll_agg{mmm} = wtAll;
    feature_mu_agg{mmm} = feature_mu;
       cluster_objects{mmm} = labels;
    
  
            
                      
                
                
    %% compute feature-and cluster space metrics
   

    % inds_labels_kmeans
    num_clusters = max(labels);
    num_markers = numel(cluster_marker_inds);
    
    %% look at feature vals
    feature_mu_reshaped = reshape(feature_mu,size(feature_mu,1),numel(fr),[]);
    % figure(19)
    % %subplot(1,2,1)
    % imagesc(1:size(agg_features,1),fr,squeeze(feature_mu_reshaped(85,:,:)))
    % set(gca,'XTick',1:size(agg_features,1),'XTickLabels',feature_labels);
    % ylabel('frequency (Hz)')
    % xlabel('Marker')
    
    % feature_sigma_reshaped = reshape(feature_sigma,size(feature_mu,1),numel(fr),...
    %     size(feature_mu,2)./numel(fr),numel(fr),size(feature_mu,2)./numel(fr));
    %
    % %subplot(1,2,2)
    % %imagesc(1:size(agg_features,1),fr,squeeze(feature_mu_reshaped(32,:,:)))
    % imagesc(squeeze(feature_sigma(32,:,:)))
    % %set(gca,'XTick',1:size(agg_features,1),'XTickLabels',feature_labels);
    % ylabel('frequency (Hz)')
    % xlabel('Marker')
    
    %% compute clustering metrics in real space
    average_pose = zeros(3,num_markers,num_markers,num_clusters);
    average_pose_struct = zeros(2,num_markers,num_markers,num_clusters);
    
    cluster_lengths = zeros(4,num_clusters);
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
    
    
    %end
 
    
    
    %number_of_subclusters = 1;
    %mmm=1;
    clustering_window = cluster_fps;
    clustering_overlap = cluster_fps./2;
    %% visualize the result of clustering
end



save_cluster_workspace=0;
if (save_cluster_workspace)
    savedirectory_subcluster =strcat(savedirectory, save_tag,filesep);
    mkdir(savedirectory_subcluster);
    save(strcat(savedirectory_subcluster ,'clusterworkspace_2.m'),'cluster_objects','clustering_ind_agg','feature_labels_agg','labels_agg','feature_mu_agg','-v7.3')
end


do_modular_plots = 0;
if (do_modular_plots)
    fprintf('starting visualization \n')
    global_cluster = 7;
    for mmm= cluster_here%1:(number_of_subclusters-1)
        figure(500+mmm)
        
        [~,sorted_ind] = sort(cluster_objects{mmm},'ASCEND');
        figure(500+mmm)
        subplot(3,1,1)
        imagesc(wtAll_agg{mmm}(sorted_ind,:)')
        caxis([0 0.1])
        set(gca,'YTick',1:numel(fr):size(wtAll_agg{mmm},1),'YTickLabels',feature_labels_agg{mmm});
        xlabel('Time (frames at 100 fps)')
        
        subplot(3,1,2)
        plot(cluster_objects{mmm}(sorted_ind))
        xlim([1 numel(cluster_objects{mmm})])
        
        subplot(3,1,3)
        
        imagesc(wtAll_agg{global_cluster}(sorted_ind,:)')
        caxis([0 0.1])
        set(gca,'YTick',1:numel(fr):size(wtAll_agg{global_cluster},1),...
            'YTickLabels',feature_labels_agg{global_cluster});
        xlabel('Time (frames at 100 fps)')
        
    end
    
    
    
    [~,sorted_ind] = sort(cluster_objects{global_cluster},'ASCEND');
    figure(500)
    subplot(3,1,1)
    imagesc(wtAll_agg{global_cluster}(sorted_ind,:)')
    caxis([0 0.1])
    set(gca,'YTick',1:numel(fr):size(wtAll_agg{global_cluster},1),'YTickLabels',feature_labels_agg{global_cluster});
    xlabel('Time (frames at 100 fps)')
    
    subplot(3,1,2)
    plot(cluster_objects{global_cluster}(sorted_ind))
    xlim([1 numel(cluster_objects{global_cluster})])
    
    agg_label_here = [];
    for mmm= 1:(number_of_subclusters-1)
        if (mmm == 1)
            agg_label_here = cluster_objects{mmm}(sorted_ind);
        else
            agg_label_here =cat(2,agg_label_here,cluster_objects{mmm}(sorted_ind));
        end
    end
    
    subplot(3,1,3)
    imagesc(agg_label_here')
    
    
    
    figure(499)
    subplot(3,1,1)
    imagesc(wtAll_agg{global_cluster}(sorted_ind,:)')
    caxis([0 0.1])
    set(gca,'YTick',1:numel(fr):size(wtAll_agg{global_cluster},1),'YTickLabels',feature_labels_agg{global_cluster});
    xlabel('Time (frames at 100 fps)')
    
    subplot(3,1,2)
    plot(cluster_objects{global_cluster}(sorted_ind))
    xlim([1 numel(cluster_objects{global_cluster})])
    
    
    % can be sped up by removing concat
    agg_label_here = [];
    agg_label_here_sorted = [];
    % subcluster_names
    tickstart = [];
    for mmm= cluster_here%1:(number_of_subclusters-1)
        tickstart = cat(1,tickstart,1+size(agg_label_here_sorted,1));
        
        %  temp = zeros(number_of_clust(mmm),numel(sorted_ind));
        temp = zeros(number_of_clust(mmm),numel(sorted_ind));
        
        %   temp(cluster_objects{mmm}(sorted_ind) ) = 1;
        
        temp(sub2ind(size(temp),cluster_objects{mmm}(sorted_ind)',1:numel(sorted_ind))) = 1;
        
        agg_label_here =cat(1,agg_label_here,temp);
        %number_of_clust(mmm)
        
        temp_sort = [];
        for jj =1:number_of_clust(global_cluster)
            globalsubind = find(cluster_objects{global_cluster}==jj);
            temp = zeros(number_of_clust(mmm),numel(globalsubind));
            [~,sorted_subind] =  sort( cluster_objects{mmm}(globalsubind),'ASCEND');
            temp(sub2ind(size(temp),...
                cluster_objects{mmm}(globalsubind(sorted_subind))',1:numel(sorted_subind))) = 1;
            temp_sort = cat(2,temp_sort,temp);
        end
        
        agg_label_here_sorted =cat(1,agg_label_here_sorted,temp_sort);
        
    end
    %conv_labels =conv(agg_label_here_sorted,ones(5,1)./5,'same');
    
    subplot(3,1,3)
    imagesc(agg_label_here_sorted)
    set(gca,'YTick',tickstart,'YTickLabels',subcluster_names(1:(number_of_subclusters-1)))
    
    figure(497)
    imagesc(agg_label_here_sorted)
    set(gca,'YTick',tickstart,'YTickLabels',subcluster_names(1:(number_of_subclusters-1)))
    
    
    cluster_to_align=2;
    feature_similarity_score = zeros(number_of_subclusters,number_of_clust(cluster_to_align));
    feature_similarity_fano = zeros(number_of_subclusters,number_of_clust(cluster_to_align));
    
    for jk = (1:number_of_subclusters)
        for ll = 1:number_of_clust(cluster_to_align)
            globalsubind = find(cluster_objects{cluster_to_align}==ll);
            cluster_mean = mean(wtAll_agg{jk}(globalsubind,:),1);
            cluster_std = std(wtAll_agg{jk}(globalsubind,:),[],1);
            feature_similarity_fano(jk,ll) = mean(cluster_std./abs(cluster_mean));
            
            feature_similarity_score(jk,ll) = mean(mean(pdist2(wtAll_agg{jk}(globalsubind(1:5:end),:),cluster_mean,...
                'euclidean')))./size(wtAll_agg{jk},2);
            
            % feature_similarity_score(jk,ll) = mean(mean(pdist(wtAll_agg{jk}(globalsubind,:),'correlation')));
        end
    end
    
    
    for jk = (1:number_of_subclusters)
        %xbins = 0:0.05:1;
        figure(600+jk)
        title(subcluster_names{jk})
        subplot(2,1,1)
        hist(  feature_similarity_score(jk,:),25)
        title(subcluster_names{jk})
        
        xlabel('Intra-cluster correlation')
        subplot(2,1,2)
        hist(  feature_similarity_fano(jk,:),25)
        title(subcluster_names{jk})
        xlim([0 2])
        
        xlabel('Intra-cluster fano')
    end
end

for mmm= 2% cluster_here%1:(number_of_subclusters-1)
    clustering_ind =   clustering_ind_agg{mmm};
    savedirectory_subcluster =strcat(savedirectory, save_tag,num2str(mmm),filesep);
    mkdir(savedirectory_subcluster);
    
    %% make spectrograms (same for all sub-marker sets)
    [~,reordered_ind] =sort(cluster_objects{mmm},'Ascend');
    feature_vis_1 = squeeze(marker_velocity(1,clustering_ind,4));
    feature_vis_2 = squeeze(marker_velocity(5,clustering_ind,4));
    feature_vis_3 = squeeze(marker_velocity(10,clustering_ind,4));
    feature_vis_4 = squeeze(markers_preproc.(marker_names{1})(clustering_ind,3));
    %
    %     [~,~,times_here,spectrogram_here] = spectrogram(  feature_vis_1(reordered_ind) ,clustering_window,clustering_overlap,cluster_fps,cluster_fps);
    %     [~,~,~,spectrogram_here_2] = spectrogram(feature_vis_2(reordered_ind),clustering_window,clustering_overlap,cluster_fps,cluster_fps);
    %     [~,~,~,spectrogram_here_3] = spectrogram(feature_vis_3(reordered_ind),clustering_window,clustering_overlap,cluster_fps,cluster_fps);
    %     [~,~,~,spectrogram_here_4] = spectrogram(feature_vis_4(reordered_ind),clustering_window,clustering_overlap,cluster_fps,cluster_fps);
    %
    %                [~,~,~,spectrogram_here_4] = spectrogram(marker_velocity(13,:,4),clustering_window,clustering_overlap,fps,fps);
    
    %  [~,~,~,spectrogram_here] = spectrogram(squeeze(delta_markers(2,3,1:frame_length_here,3)),floor(fps),floor(fps./2),fps,fps);
    time_ordering = [];
    cluster_vals = [];
    time_ordering_fulltrace = cell(1,num_clusters);
    cluster_ordering_fulltrace = cell(1,num_clusters);
    
    pose_ordering = cell(numel(marker_names),numel(marker_names));
    pose_ordering_full = cell(numel(marker_names),numel(marker_names));
    time_axis_full =  zeros(1,numframes );
    cluster_size =zeros(1,num_clusters);
    
    %% initialize
    
    for jj = 1:numel(marker_names)
        for lk = (jj+1):num_markers
            pose_ordering_full{jj,lk} = zeros(1,numframes );
            %  pose_ordering{jj,lk}
        end
    end
    head_height = zeros(1,numframes );
    
    %% get positions in xyz for markers in each cluster
    fprintf('getting cluster temporal positions \n')
    lastind = zeros(num_markers,num_markers);
    lastclust = zeros(num_markers,num_markers);
    cluster_quality = zeros(3,num_clusters);
    cluster_size = zeros(1,num_clusters);
    
    
    
    for ll = 1:number_of_clust(mmm)
        %% ordering of the spectrogram
        fprintf('cluster %f \n',ll)
        time_ordering = cat(1,time_ordering,find(cluster_objects{mmm}==ll));
        cluster_vals = cat(1,cluster_vals,ll*ones(numel(find(cluster_objects{mmm}==ll)),1));
        cluster_size(ll) = numel(find(cluster_objects{mmm}==ll));
        
        time_ordering_fulltrace{ll} = find(cluster_objects{mmm}==ll);
        %             (unique(bsxfun(@plus,round(times_here(find(cluster_objects{mmm}==ll))*clustering_window)',...
        %                 -floor(clustering_overlap):floor(clustering_overlap))));
        
        %  cluster_ordering_fulltrace{ll} = repmat(
        
        % (bsxfun(@plus,round(times_here(find(cluster_objects{mmm}==ll))*clustering_window)',...
        %    -floor(clustering_overlap):floor(clustering_overlap)));
        
        time_ordering_fulltrace{ll}(time_ordering_fulltrace{ll}<1) = 1;
        
        %time_ordering_fulltrace{ll} = time_ordering_fulltrace{ll}+min( clustering_ind);
        
        time_ordering_fulltrace{ll}(time_ordering_fulltrace{ll} > max( clustering_ind)) = max( clustering_ind);
        time_ordering_fulltrace{ll}=   clustering_ind(time_ordering_fulltrace{ll});
    end
    
    
    good_clusters = find(cluster_size>=100);%intersect(find(cluster_size<25),find(cluster_size>=10)); %,
    movie_output = cell(1,numel(good_clusters));
    
    frames_use_anim = 300;
    movies_to_examine = good_clusters;% [1,9,17,18,19,39,41,43,95,89,50,51,59];%good_clusters(good_clusters>47);%[83,60,3,6,28,35,59,42,71,83];%[3,6,23];
    save_movie = 1;
    %save_tags = {'1_torque','9_raiseexplore','17extend','18drop1','19rise','39walknshake','41drop2','43uppershake','95walkshake','89slowwalk','50rearinvestigate','51slowshift','59turn'};
    %save_tags = {'KX83_ketamine_wobble_walk','KX60_FullTiltWobble','KX3_wobblerear','KX6_wobblewalk_headmove',...
    %    'KX28_slow_odor_search','KX35_TurnLeft','KX59_turnright','KX42wobbleup','KX71_tighttwist','KX83_fastodorsearch'};
    
    save_tags_temp = (num2str(good_clusters'));
    for kk =1:numel(good_clusters)
        save_tags{kk} = save_tags_temp(kk,:);
    end
    
    
    plot_spectrogram =0;
    if plot_spectrogram
        for jjj = good_clusters(1:end);%58;%36;%good_clusters(1:end)
            figure(470)
            % title(strcat('cluster ',num2str(jjj)))
            %      frame_inds = time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use,numel(time_ordering_fulltrace{jjj})));
            %          temp = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
            %  M = [];
            %   movie_output{jjj} = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
            
            trace_1 = bsxfun(@minus,squeeze(delta_markers(3,6,(time_ordering_fulltrace{jjj}),4)),...
                mean(squeeze(delta_markers(3,6,(time_ordering_fulltrace{jjj}),4))));
            
            trace_2 = bsxfun(@minus,squeeze(delta_markers(3,8,(time_ordering_fulltrace{jjj}),4)),...
                mean(squeeze(delta_markers(3,8,(time_ordering_fulltrace{jjj}),4))));
            
            trace_3 = bsxfun(@minus,squeeze(delta_markers(4,5,(time_ordering_fulltrace{jjj}),4)),...
                mean(squeeze(delta_markers(4,5,(time_ordering_fulltrace{jjj}),4))));
            
            trace_4 = bsxfun(@minus,squeeze(delta_markers(3,6,(time_ordering_fulltrace{jjj}),4)),...
                mean(squeeze(delta_markers(3,6,(time_ordering_fulltrace{jjj}),4))));
            
            [~,~,times_here,spectrogram_here_clust] = spectrogram(trace_1,...
                clustering_window,clustering_overlap,fps,fps);
            [~,~,~,spectrogram_here_2_clust] = spectrogram(trace_2,...
                clustering_window,clustering_overlap,fps,fps);
            [~,~,~,spectrogram_here_3_clust] = spectrogram(trace_3,...
                clustering_window,clustering_overlap,fps,fps);
            
            [~,~,times_here,spectrogram_here_4_clust] = spectrogram(squeeze(marker_velocity(1,(time_ordering_fulltrace{jjj}),4)),...
                clustering_window,clustering_overlap,fps,fps);
            
            color_axis = [-8 0];
            freq_limits = [0 40];
            %xlimits = [1 500];
            % xlimits = [1 numel(cluster_vals)];
            uper_freq = 75;
            subplot(4,2,1)
            imagesc(log10(spectrogram_here_clust))
            title('head to spine')
            caxis( color_axis)
            ylim( freq_limits)
            % xlim(xlimits)
            ylabel('Frequency [Hz]')
            
            subplot(4,2,3)
            imagesc(log10(spectrogram_here_2_clust))
            title('hip to spine')
            caxis( color_axis)
            ylim( freq_limits)
            % ylabel('Frequency [Hz]')
            xlabel('time ')
            % xlim(xlimits)
            
            subplot(4,2,5)
            imagesc(log10(spectrogram_here_3_clust))
            title('inter spine')
            caxis( color_axis)
            ylim( freq_limits)
            % xlim(xlimits)
            %    ylabel('Frequency [Hz]')
            
            subplot(4,2,7)
            imagesc(log10(spectrogram_here_4_clust))
            title('head velocity')
            caxis( color_axis)
            ylim( freq_limits)
            % xlim(xlimits)
            %ylabel('Frequency [Hz]')
            
            
            %  figure(366)
            axis_ordered = 0:1./fps:(numel(pose_ordering_full{2,3})./fps-1./fps);
            
            subplot(4,2,2)
            plot(squeeze(delta_markers(3,4,time_ordering_fulltrace{jjj},4)))
            title('head to spine')
            ylabel('marker distance (mm)')
            
            subplot(4,2,4)
            plot(squeeze(delta_markers(3,8,time_ordering_fulltrace{jjj},4)))
            title('head to hip')
            %     ylabel('marker distance (mm)')
            
            subplot(4,2,6)
            plot(squeeze(delta_markers(4,5,time_ordering_fulltrace{jjj},4)))
            title('inter spine')
            %  ylabel('marker distance (mm)')
            
            subplot(4,2,8)
            plot(245*squeeze(marker_velocity(1,time_ordering_fulltrace{jjj},4)))
            title('marker velocity')
            ylabel('head velocity (mm/s)')
            ylim([0 250])
            
            print('-depsc',strcat(savedirectory_subcluster,'plotsfor',num2str(jjj),'.eps'))
            print('-dpng',strcat(savedirectory_subcluster,'plotsfor',num2str(jjj),'.png'))
            
            %             figure(470)
            %                         plot(pose_ordering_full{2,3}...
            %                             (bsxfun(@plus,round(times_here(find(cluster_objects{mmm}==jjj))*clustering_window)',...
            %                 -floor(clustering_overlap):floor(clustering_overlap)))')
            
            %             ylim([0 100])
            %
            %             (bsxfun(@plus,round(times_here(find(cluster_objects{mmm}==ll))*clustering_window)',...
            %                 -floor(clustering_overlap):floor(clustering_overlap)))
            %
        end
    end
    
    
    
    %
    %                if (save_movie)
    %
    %
    %             end
    feature_mu_reshape = reshape(feature_mu_agg{mmm},size(feature_mu_agg{mmm},1),numel(fr),[]);
    
    x_ind = size(feature_mu_reshape,2);
    y_ind = size(feature_mu_reshape,3);
    
    xcorr_features = zeros(numel(good_clusters),numel(good_clusters));
    euc_features = zeros(numel(good_clusters),numel(good_clusters));
    wave_features = zeros(numel(good_clusters),numel(good_clusters));
    mean_wave = zeros(max(good_clusters),size(wtAll_agg{mmm},2));
    for lkl = 1:max(good_clusters)
        mean_wave(lkl,:) = mean(wtAll_agg{mmm}(find(cluster_objects{mmm}==lkl),:),1);
    end
    
    
    for jjj = good_clusters(1:end);
        for mkm = good_clusters(find(good_clusters == jjj):end)
            corr_c = xcorr2(squeeze(feature_mu_reshape(jjj,:,:)),squeeze(feature_mu_reshape(mkm,:,:)));
            xcorr_features(jjj,mkm) = corr_c(x_ind,y_ind);
            euc_features(jjj,mkm) = pdist2(feature_mu_agg{mmm}(jjj,:),feature_mu_agg{mmm}(mkm,:),'correlation');
            
            wave_features(jjj,mkm)  = pdist2(mean_wave(jjj,:),...
                mean_wave(mkm,:),'correlation');
        end
    end
    xcorr_t = xcorr_features';
    xcorr_features(find(tril(ones(size( xcorr_features)),0))) = xcorr_t(find(tril(ones(size( xcorr_features)),0)));
    
    euc_t = euc_features';
    euc_features(find(tril(ones(size( euc_features)),0))) = euc_t(find(tril(ones(size(euc_features)),0)));
    
    wave_t = wave_features';
    wave_features(find(tril(ones(size( wave_features)),0))) = wave_t(find(tril(ones(size( xcorr_features)),0)));
    %
    [idx,C,sumd,D] = kmeans(feature_mu_agg{mmm}(good_clusters,:),12,'distance','correlation');
    % [idx,C,sumd,D] = kmeans( mean_wave(good_clusters,:),12,'distance','correlation');
    
    [vals,inds] = sort(idx,'ASCEND');
    
    
    figure(371)
    imagesc(xcorr_features(good_clusters(inds),good_clusters(inds)))
    caxis([-10 10])
    print('-depsc',strcat(savedirectory_subcluster,'xcorrclustered',num2str(mmm),'.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'xcorrclustered',num2str(mmm),'.png'))
    
    figure(373)
    imagesc(euc_features(good_clusters(inds),good_clusters(inds)))
    caxis([0 3])
    print('-depsc',strcat(savedirectory_subcluster,'eucclustered',num2str(mmm),'.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'eucclustered',num2str(mmm),'.png'))
    
    figure(372)
    imagesc(wave_features(good_clusters(inds),good_clusters(inds)))
    caxis([0 0.4])
    
    
    
    
    good_clusters_sorted=good_clusters(inds);
    
    
    
    do_cluster_plots = 1;
    if (do_cluster_plots)
        
        for jjj =  good_clusters_sorted(1:end)%1:end);%58;%36;%good_clusters(1:end)
            
            
            save_ind = find(good_clusters_sorted(1:end)== jjj);
            
            figure(380)
            %subplot(1,2,1)
            
            imagesc(1:size(feature_mu_reshape,3),fr,squeeze(feature_mu_reshape(jjj,:,:)))
            set(gca,'XTick',1:3:3*size(feature_mu_reshape,3),'XTickLabels',feature_labels_agg{mmm});
            ylabel('frequency (Hz)')
            xlabel('Marker')
            title('Feature weights')
            print('-depsc',strcat(savedirectory_subcluster,'coeffweightsfor',num2str(save_ind),'.eps'))
            print('-dpng',strcat(savedirectory_subcluster,'coeffweightsfor',num2str(save_ind),'.png'))
            
            
            %
            % figure(400)
            % %imagesc(squeeze(average_pose(2,:,:,jjj))./squeeze(average_pose(1,:,:,jjj)))
            % set(gca,'XTick',1:size(average_pose,2),'XTickLabels',marker_names);
            % set(gca,'YTick',1:size(average_pose,2),'YTickLabels',marker_names);
            % title('intermarker std./global std')
            
            
            
            figure(390)
            imagesc(wtAll_agg{mmm}(find(cluster_objects{mmm}==jjj),:)')
            caxis([0 2])
            set(gca,'YTick',1:3*numel(fr):size(wtAll_agg{mmm},1),'YTickLabels',feature_labels_agg{mmm});
            xlabel('Time (frames at 100 fps)')
            title('Observed wavelet coeffs')
            
            print('-depsc',strcat(savedirectory_subcluster,'coeffplotsfor',num2str(save_ind),'.eps'))
            print('-dpng',strcat(savedirectory_subcluster,'coeffplotsfor',num2str(save_ind),'.png'))
        end
        
        
        
        
        
        dH = designfilt('highpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', 1/(fps/2), ...
            'DesignMethod', 'butter');
        [f1,f2] = tf(dH);
        abs_velocity_antialiased_hipass = abs_velocity_antialiased;
        %    for mm = 1:size(abs_velocity_antialiased_hipass,1)
        %        abs_velocity_antialiased_hipass(mm,:)  = filtfilt(f1,f2,squeeze(marker_position(mm,:,3)));
        %    end
        
        for jjj =good_clusters_sorted(1:end)%1:end);%58;%36;%good_clusters(1:end)
            
            
            save_ind = find(good_clusters_sorted(1:end)== jjj);
            
            %label individual instances
            inst_label = zeros(1,numel(clustering_ind));
            inst_label(cluster_objects{mmm} == jjj) = 1;
            pixellist = bwconncomp(inst_label);
            
            num_rand = 40;
            marker_align = 3;%cluster_markersets{mmm}(1);
            markers_to_plot = [3]%,5*3,find(cluster_markersets{mmm} == 17)*3];
            amt_extend_xcorr = 75; %in 3x downsample
            amt_extend_plot = 80;
            num_clust_sample = min(pixellist.NumObjects,  num_rand );
            
            if num_clust_sample
                rand_inds = randsample(1:pixellist.NumObjects,...
                    num_clust_sample);
                
                clust_ind_here = cell(1,num_clust_sample);
                trace_length_here = 2*amt_extend_xcorr*3+1;
                
                av_trace = zeros(1,trace_length_here);
                
                for lk =1:num_clust_sample
                    % clustmed = median(pixellist.PixelIdxList{rand_inds(lk)})
                    indexh = unique(bsxfun(@plus,floor(median(pixellist.PixelIdxList{rand_inds(lk)})),-amt_extend_xcorr:amt_extend_xcorr));
                    indexh(indexh<1) = 1;
                    indexh(  indexh>numel(clustering_ind)) = numel(clustering_ind);
                    indexh = unique(indexh);
                    clusterind = unique(clustering_ind(indexh));
                    clust_ind_here{lk}  =clusterind;
                    % abs_velocity_antialiased
                    trace_add = agg_features(marker_align,(clusterind(1):clusterind(end)) )./num_clust_sample;
                    
                    
                    av_trace(1:min(length(trace_add), trace_length_here)) = av_trace(1:min(length(trace_add), trace_length_here))...
                        +  trace_add(1:min(length(trace_add), trace_length_here))./num_clust_sample;
                end
                %  av_trace = av_trace+marker_velocity(marker_align,clusterind(1):clusterind(end) ,4)./num_clust_sample;
                
                
                for llll=1:2
                    av_trace = zeros(1,trace_length_here);
                    for lk =1:num_clust_sample
                        %[c,lags] =   xcorr(marker_velocity(marker_align,clust_ind_here{lk}(1):clust_ind_here{lk}(end) ,4),av_trace,amt_extend_xcorr);
                        [c,lags] =   xcorr(agg_features(marker_align,clust_ind_here{lk}(1):clust_ind_here{lk}(end) ),av_trace,amt_extend_xcorr);
                        
                        [~,ind] = max(c);
                        %   lags(ind) = 0;
                        ind
                        clust_ind_here{lk} = bsxfun(@plus,clust_ind_here{lk}, lags(ind));
                        % av_trace = av_trace+marker_velocity(marker_align,clust_ind_here{lk}(1):clust_ind_here{lk}(end) ,4);
                        av_trace(1:min(length(trace_add), trace_length_here)) = av_trace(1:min(length(trace_add), trace_length_here))+...
                            trace_add(1:min(length(trace_add), trace_length_here))./num_clust_sample;
                        
                    end
                end
                
                
                y_gap = 3;
                for mm = 1:numel(markers_to_plot)
                    figure(500+mm)
                    velocity_size = size(abs_velocity_antialiased_hipass,2);
                    for lk =1:min(pixellist.NumObjects,num_rand)
                        % subplot(2,numel(markers_to_plot)./2,mm)
                        %   subplot(1,1,mm)
                        velocity_plot = y_gap*lk+agg_features(   markers_to_plot(mm), max(1,(clust_ind_here{lk}(1)-amt_extend_plot)):...
                            min(velocity_size,clust_ind_here{lk}(end)+amt_extend_plot));
                        
                        %                        position_plot = marker_position(   markers_to_plot(mm), max(1,(clust_ind_here{lk}(1)-amt_extend_plot)):...
                        %                            min(velocity_size,clust_ind_here{lk}(end)+amt_extend_plot),3);
                        %                        position_plot = bsxfun(@minus,position_plot,mean(position_plot));
                        %
                        %    plot(marker_velocity(   markers_to_plot(mm), (clust_ind_here{lk}(1)-amt_extend_plot):(clust_ind_here{lk}(end)+amt_extend_plot)  ,4))
                        plot(velocity_plot )
                        
                        
                        hold on
                    end
                    title(marker_names{floor(markers_to_plot(mm)./3)})
                    box off
                    max_time_val = (clust_ind_here{lk}(end)+amt_extend_plot)-(clust_ind_here{lk}(1)-amt_extend_plot);
                    xtickshere = 1:25:max_time_val;
                    xlim([1 max_time_val])
                    set(gca,'XTick',xtickshere,'XTickLabel',num2str(floor(1000*xtickshere./300)'))
                    ylabel('Marker velocity (mm/frame)')
                    xlabel('Time (ms)')
                    Lhere =  get(gca,'ylim');
                    %    ylim([0 20])
                    ylim([-20 min(Lhere(2),120)])
                    hold off
                end
                print(strcat(savedirectory_subcluster,'timetraces',num2str(save_ind),'.pdf'),'-dpdf','-bestfit')
                print('-dpng',strcat(savedirectory_subcluster,'timetraces',num2str(save_ind),'.png'))
                
                %
                %                 figure(480)
                %             plot(markers_preproc.HeadF(clust_ind_here{lk}(1):clust_ind_here{lk}(end),3))
                %              hold on
                % %
                
            end
        end
    end
    
    
    clipped_clustering_ind = clipped_index_agg{mmm}(clustering_ind);
    do_movies =1;
    num_ex = 10;
    if (do_movies)
        for jjj =  good_clusters_sorted(3:end)%1:end);%58;%36;%good_clusters(1:end)
            
            
            save_ind = find(good_clusters_sorted(1:end)== jjj);
            
            matlab_fr = 2;
            amt_extend = 50;
            
            figure(370)
            align_clusters = 1;
            if align_clusters
                %  movie_output{jjj} = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
                
                
                %label individual instances
                inst_label = zeros(1,numel(clipped_clustering_ind));
                inst_label(cluster_objects{mmm} == jjj) = 1;
                pixellist = bwconncomp(inst_label);
                
                indexh = unique(bsxfun(@plus,find(cluster_objects{mmm} == jjj),- amt_extend: amt_extend));
                indexh(indexh<1) = 1;
                indexh(  indexh>numel(clipped_clustering_ind)) = numel(clipped_clustering_ind);
                indexh= unique(indexh);
                
                clusterind_full =  unique(clipped_clustering_ind(indexh));
                cluster_mean =  nanmean(markers_preproc.SpineM( clusterind_full ,:),1);
                cluster_mean_array =  (markers_preproc.SpineM( clusterind_full ,:));
                
                
                
                markers_preproc_cluster_aligned = struct();
                markernames = fieldnames(markers_preproc);
                
                fprintf('rotating markers')
                rotangle = atan2(-(markers_preproc.SpineF(clusterind_full,2)-markers_preproc.SpineM(clusterind_full,2)),...
                    (markers_preproc.SpineF(clusterind_full,1)-markers_preproc.SpineM(clusterind_full,1)));
                global_rotmatrix = zeros(2,2,numel(rotangle));
                global_rotmatrix(1,1,:) = cos(rotangle);
                global_rotmatrix(1,2,:) = -sin(rotangle);
                global_rotmatrix(2,1,:) = sin(rotangle);
                global_rotmatrix(2,2,:) = cos(rotangle);
                
                M = [];
                
                %look at at most 50 clusters
                
                
                figure(370)
                
                %  v = VideoWriter(strcat(savedirectory_subcluster,'movie',save_tags{find(movies_to_examine==jjj)}),'MPEG-4');
                v = VideoWriter(strcat(savedirectory_subcluster,'movie',num2str(save_ind)),'MPEG-4');
                
                open(v)
                
                rand_inds = randsample(1:pixellist.NumObjects,...
                    min(pixellist.NumObjects,num_ex));
                
                for lk =1:min(pixellist.NumObjects,num_ex)
                    
                    indexh = unique(bsxfun(@plus,pixellist.PixelIdxList{rand_inds(lk)},-amt_extend:amt_extend));
                    indexh(indexh<1) = 1;
                    indexh(  indexh>numel(clipped_clustering_ind)) = numel(clipped_clustering_ind);
                    indexh = unique(indexh);
                    
                    clusterind = unique(clipped_clustering_ind(indexh));
                    
                    clustersubind = arrayfun(@(x)(find(clusterind_full == x)),clusterind);
                    
                    for mm = 1:numel(markernames);
                        markers_preproc_cluster_aligned.(markernames{mm}) = markers_preproc.(markernames{mm})(clusterind_full(clustersubind),:);
                        markers_preproc_orig.(markernames{mm}) = markers_preproc.(markernames{mm})(clusterind_full(clustersubind),:);
                    end
                    
                    
                    for mm = (1:numel(markernames));
                        markers_preproc_cluster_aligned.(markernames{mm}) =  bsxfun(@minus,...
                            markers_preproc_cluster_aligned.(markernames{mm}), cluster_mean_array(clustersubind,:));
                    end
                    
                    %% now make and apply a rotation matrix
                    
                    
                    for mm = (1:numel(markernames));
                        markers_preproc_cluster_aligned.(markernames{mm})(:,1:2) =  ...
                            squeeze(mtimesx(global_rotmatrix(:,:,clustersubind),(reshape(markers_preproc_cluster_aligned.(markernames{mm})(:,1:2)',2,1,numel(clusterind)))...
                            ))';
                        
                        markers_preproc_cluster_aligned.(markernames{mm}) = bsxfun(@plus,markers_preproc_cluster_aligned.(markernames{mm}), cluster_mean);
                    end
                    
                    
                    
                    %                 markers_preproc_anim = struct();
                    %         markernames = fieldnames(markers_preproc);
                    %         for mm = 1:numel(markernames);
                    %            markers_preproc_anim.(markernames{mm}) = markers_preproc_cluster_aligned.(markernames{mm})(clusterind,:);
                    %         end
                    frame_inds = 1:matlab_fr:min(frames_use_anim,numel( clusterind));
                    h=  figure(370)
                    movie_output_temp = animate_markers_clusteraligned(markers_preproc_cluster_aligned,markers_preproc_orig,frame_inds',marker_names,markercolor,links,M,jjj,lk,h, subcluster_names{mmm});
                    
                    % writeVideo(v,movie_output{jjj})
                    writeVideo(v,movie_output_temp)
                    
                end
                close(v)
                
            else
                
                % title(strcat('cluster ',num2str(jjj)))
                frame_inds = (time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use_anim,numel(time_ordering_fulltrace{jjj}))))';
                
                %          temp = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
                M = [];
                
                
                movie_output_temp = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,M,jjj);
                
                
                v = VideoWriter(strcat(savedirectory_subcluster,'movie',save_tags{find(movies_to_examine==jjj)}),'MPEG-4');
                open(v)
                % writeVideo(v,movie_output{jjj})
                writeVideo(v,movie_output_temp)
                
                close(v)
            end
        end
    end
end


%             figure(370+mmm)
%             title(num2str(jjj))
%
%
%             plot3( squeeze(markers_preproc.(marker_names{1})(1,1)),...
%                 squeeze(markers_preproc.(marker_names{1})(1,2)),...
%                 squeeze(markers_preproc.(marker_names{1})(1,3)),'o','Color',markercolor{jj})
%             ax = gca;
%             axis(ax,'manual')
%             zlim([0 250])
%             xlim([-200 250])
%             ylim([-200 250])
%
%             matlab_fr = 10;
%
%             frame_inds = time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use,numel(time_ordering_fulltrace{jjj})));
%
%             %M{jjj} = movie;
%             %Mhere = movie;
%             for lk = frame_inds'%1:10:10000
%                 ind_to_plot = lk;
%
%
%                 set(gca,'Nextplot','ReplaceChildren');
%                 handles_here = cell(1,numel(marker_names));
%                 for jj = 1:numel(marker_names)
%                     handles_here{jj} = plot3( squeeze(markers_preproc.(marker_names{jj})(ind_to_plot,1)),...
%                         squeeze(markers_preproc.(marker_names{jj})(ind_to_plot,2)),...
%                         squeeze(markers_preproc.(marker_names{jj})(ind_to_plot,3)),'o','Color',markercolor{jj},'MarkerFaceColor','auto');
%                     hold on
%                 end
%                 for mm = 1:numel(links)
%                     plot3( [squeeze(markers_preproc.(marker_names{links{mm}(1)})(ind_to_plot,1)) ...
%                         squeeze(markers_preproc.(marker_names{links{mm}(2)})(ind_to_plot,1)) ],...
%                         [ squeeze(markers_preproc.(marker_names{links{mm}(1)})(ind_to_plot,2))...
%                         squeeze(markers_preproc.(marker_names{links{mm}(2)})(ind_to_plot,2))],...
%                         [squeeze(markers_preproc.(marker_names{links{mm}(1)})(ind_to_plot,3))...
%                         squeeze(markers_preproc.(marker_names{links{mm}(2)})(ind_to_plot,3))],'Color',markercolor{links{mm}(1)});
%                 end
%                 %delete
%                 if (save_movie)
%                     M{jjj}(find(frame_inds == lk)) = getframe(gcf);
%                 end
%
%                 drawnow
%                 hold off
%             end
%             set(gca,'Nextplot','add');




%             for jjj = movies_to_examine;%58;%36;%good_clusters(1:end)
%
%                  v = VideoWriter(strcat(savedirectory_subcluster,'slow',save_tags{find(movies_to_examine==jjj)}),'MPEG-4');
%                 open(v)
%                 writeVideo(v,movie_output{jjj})
%                 close(v)
%             end

plot_simul = 0;
if (plot_simul)
    nrows = 2;
    ncols = 3;
    
    %  cluster_numbers = [2,42,51,55,59,124,73,69,65,62,183,113,112,164,158,175,15,183];
    %             cluster_numbers = [6,50,11,43,35,33,29,22,16,13,12,10];
    cluster_numbers = [2,6,25,83,81,100];
    
    %      good_clusters(1:16);
    %            cluster_names = {'tap/drop','groom hands','odor sample','up/down/up',....
    %                'under and up','tap and rise','sniff and explore','left tap','lick',...
    %                'tap/drop','hi sample','tap','lick and up','med sample','lever sample','low sample'};
    %                     cluster_names = {'hi sample','sniff/explore','mid sniff','mid sample','low sample','head crane','reach hisample','low sniff',...
    %                         'eating','sniff low','down and up','sniff and shake','sniff and explore'};
    cluster_names = {'CCW turn','CW turn','straight','tight turn','CCW hi tilt','CW hi tilt'};
    
    
    %num2str(cluster_numbers');
    fighand =   figure(388)
    set(fighand,'Color','k')
    set(fighand,'Position',[100 100 1100 1100])
    for lk = 1:1000 %numel(movie_output{jjj})
        for ll = 1:nrows*ncols  %numel(cluster_numbers)
            mov_ind = cluster_numbers(ll);
            subplot_tight(nrows,ncols, ll)
            movie_size = numel(movie_output{mov_ind});
            frame_use = mod(lk,movie_size);
            if (frame_use==0)
                frame_use = 1;
            end
            
            imshow(movie_output{mov_ind}(frame_use).cdata,movie_output{mov_ind}(frame_use).colormap)
            title(cluster_names{ll},'Color','w')
        end
        M_here(lk) =      getframe(gcf);
        
    end
    
    v = VideoWriter(strcat(savedirectory_subcluster,'aggregate_movie_2'),'MPEG-4');
    open(v)
    writeVideo(v, M_here)
    close(v)
    clear M_here
end



%% randomly remove 1% of the data for a single marker and impute repeatedly
number_of_repeats = 10;
removal_fraction = 1; %percent
buffer_length = 20;

% first select frames not within 20 frames of a true bad frame
bad_times_buffer = bsxfun(@plus,bad_times,-buffer_length:1:buffer_length);
bad_times_buffer = unique(reshape(bad_times_buffer,1,[]));
bad_times_buffer(bad_times_buffer<1) = 1;
bad_times_buffer(bad_times_buffer>marker_frame_length) = marker_frame_length;
good_times_buffered = setxor(1:marker_frame_length,bad_times_buffer);

% find times with a certain number of frames beforehand
marker_remove = 1;
random_times = randsample(good_times_buffered,floor(0.01*removal_fraction*marker_frame_length));
%build the data table for the regression
% svmtable = table();
% for j =1:numel(marker_names)
%     for lk = 1:3
%     svmtable.(strcat(marker_names{j},num2str(lk))) = markers_preproc.(marker_names{j})(:,lk);
%     end
% end
%
%
%
% Mdl = fitrsvm(svmtable(good_times_buffered,:),strcat(marker_names{marker_remove},num2str(1)));

%% plot the data
%interesting things to plot: DISTANCES spine length (3,5) shows 4 hz oscillatory
%features, also higher freq
% time series of offset to spine mid (3,4) shows 'shake' instances
% velocity of marker 4 also shows the shake

x_axis = 0:1/fps:(marker_frame_length-1)./fps;
x_axis_analog = 0:1/analog_fps:(numel(lever_thresholded)-1)./analog_fps;


scale_factor = 25;
%plot
figure(24)
%plot(x_axis,squeeze(delta_markers(2,3,:,4)),'b')
%plot(x_axis,squeeze(delta_markers(2,3,:,4)))
reduce_plot(x_axis,markers_preproc.LHand(:,3),'r')
hold on

reduce_plot(x_axis,markers_preproc.RHand(:,3),'b')
reduce_plot(x_axis,markers_preproc.FHead(:,3),'c')

%plot(x_axis,markers_preproc.(marker_names{2})(:,3),'r')

reduce_plot(x_axis_analog,scale_factor*analog.LICK_SENSOR,'g')
reduce_plot(x_axis_analog,scale_factor*lever_thresholded,'k')
reduce_plot(x_axis_analog,scale_factor*analog.SPEAKER,'y')
ylim([-50 250])
hold off

figure(30)
plot(x_axis,squeeze(delta_markers(7,8,:,4)),'r')
hold on
plot(x_axis,squeeze(delta_markers(2,3,:,4)),'b')
hold off

matrix_here = squeeze(delta_markers(3,5,:,3));
matrix_here = squeeze(marker_velocity(4,:,4));
%matrix_here = squeeze(coeff(:,8)'*agg_features);

%spectrogram with a 1s window (fps) 50% overlap (fps/2) with 1 Hz sampling
%(Nfft = fps) and the given framerate
[S,F,T,P] = spectrogram(matrix_here,fps,floor(fps./2),fps,fps);
%
% params.tapers = [2 10];
% params.Fs = 245;
% [S,t,f]= mtspecgramc(matrix_here,[1 1],params);
figure(25)
imagesc(log10(P))

figure(26)
plot(x_axis,marker_velocity(4,:,4))


figure(23)
subplot(2,1,1)
plot(squeeze(markers.R_Head))
hold on
plot(squeeze(markers.Tail))

%plot(bsxfun(@times,(resample_analog'),[1000 10^4]')')
hold off
axis tight
subplot(2,1,2)
%plot(squeeze(analog.LEVER))
plot(squeeze(markers.L_Foot))
hold on
plot(squeeze(markers.R_Foot))

%plot(squeeze(analog.SPEAKER))

axis tight

time_inds = 1:5000;
figure(24)
plot3(markers.R_Head(time_inds ,1),markers.R_Head(time_inds ,2),markers.R_Head(time_inds ,3))
hold on
plot3(markers.Tail(time_inds ,1),markers.Tail(time_inds ,2),markers.Tail(time_inds ,3),'r')
hold off

find(analog.SPEAKER==0,1,'first')