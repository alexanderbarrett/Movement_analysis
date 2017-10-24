%% GMM clustering on the spectrograms

function [cluster_struct] = Cluster_GMM_clean(agg_features,opts,frames)



%% set the number of voices/octaves for the GMM
numVoices = 10;
numOctave = 6;

%% set wavelet properties
wname = 'gaus1';
dt = 0.01;
f0 = centfrq(wname);

%% set the frequency for the scales -- these are used later on tog enerate the 
minfreq = 0.5;
minscale = f0./(minfreq*dt);
s0 = minscale*dt;

a0 = 2^(1/numVoices);
scales = s0*a0.^(0:numOctave*numVoices);

%% get the frequencies from the wavelet -- this is what you ultimately want
Freq = scal2frq(scales,wname,dt);

pc_spectrograms = cell(1,size(agg_features,1));

for k=1:size(agg_features,1)
    
   %% get a multiresolution spectrogram
    fprintf('starting multiresolution spectrogram for feature %f \n',k);
    fr = [];
    freq_rev = fliplr(Freq);
    for mm =1:numOctave
        freq_subset = freq_rev(1+(mm-1)*numVoices:(mm)*numVoices);
        freq_delta_here = (max(freq_subset)-min(freq_subset))./numVoices;
        freq_range_here = min(freq_subset):freq_delta_here:max(freq_subset);
        %      freq_range_here
        [~,fr_temp,time_clustering,pc_spectrograms_temp] = spectrogram(agg_features(k,:),opts.clustering_window,...
            opts.clustering_overlap,freq_range_here,opts.fps);
        pc_spectrograms{k} =   cat(1,pc_spectrograms{k},pc_spectrograms_temp);
        fr = cat(1,fr,fr_temp);
    end
   
end
num_fr = numel(fr);
agg_spectrograms = cell2mat(pc_spectrograms'); %second dimension is time base

%% normalize the spectrograms
agg_spectrograms = log10(agg_spectrograms);
agg_spectrograms(isinf(agg_spectrograms)) = -20;
agg_spectrograms(isnan(agg_spectrograms)) = -20;

timebase = size(pc_spectrograms{1},2);

frames_per_bin = floor(size(agg_features,2)./timebase);

%% PCA on spectrograms
[coeff2,score2,latent2,tsquared2,explained2] = pca(agg_spectrograms(:,:)');
pca_agg_spectrogram =(score2');


fprintf('number of PCS2 %f size of agg spectrogram %f \n',opts.num_pcs_2,size(pca_agg_spectrogram,1));
opts.num_pcs_2 = min(size(pca_agg_spectrogram,1),opts.num_pcs_2);

%% do a GMM clustering
obj = fitgmdist((pca_agg_spectrogram(1:opts.num_pcs_2 ,:))',opts.num_clusters,...
    'Start', 'plus', 'Options', statset('MaxIter',1000),'RegularizationValue',0.1);
clusterobj = cluster(obj,(pca_agg_spectrogram(1:opts.num_pcs_2 ,:))');

labels =clusterobj;
featurespectrogram = pca_agg_spectrogram(1:opts.num_pcs_2 ,:);


feat_pcs_full = reshape(coeff2,num_fr,size(agg_features,1),size(coeff2,2));
cluster_struct.feat_pcs = feat_pcs_full(:,:,1:opts.num_pcs_2);

V = coeff2(:,1:opts.num_pcs_2);

mu = obj.mu;
sigma = obj.Sigma;

feature_mu = mu * V';
feature_sigma = zeros(size(coeff2,1),size(coeff2,1) ,opts.num_clusters);
for j = 1:opts.num_clusters
    feature_sigma(:,:,j) = V * sigma(:,:,j) * V';
end


%% convert the labels to be on each frame

labels_longer = zeros(1,numel(frames));
wtAll = zeros(size(featurespectrogram,1),numel(frames));

for jjj = 1:max(labels)
    ind_label = unique(reshape(bsxfun(@plus,round(time_clustering(find(labels==jjj))*opts.clustering_window)',...
        -floor(opts.clustering_overlap):floor(opts.clustering_overlap)),1,[]));
    ind_label(ind_label <1) = 1;
    ind_label(ind_label>numel(frames)) = numel(frames);
    labels_longer(ind_label) = jjj;
end

for jjj = 1:size(featurespectrogram,2)
    ind_label = unique(reshape(bsxfun(@plus,round(time_clustering(jjj)*opts.clustering_window)',...
        -floor(opts.clustering_overlap):floor(opts.clustering_overlap)),1,[]));
    ind_label(ind_label <1) = 1;
    ind_label(ind_label>numel(frames)) = numel(frames);
    wtAll(:,ind_label) = repmat(featurespectrogram(:,jjj),1,numel(ind_label));
end

%% output a cluster structure
cluster_struct = struct();
cluster_struct.labels = labels_longer;
cluster_struct.labels_orig = labels;
cluster_struct.num_clusters = length(unique(cluster_struct.labels));
cluster_struct.num_clusters_req = opts.num_clusters;
cluster_struct.wtAll = wtAll';
cluster_struct.fr = fr;

cluster_struct.feature_mu = feature_mu;
cluster_struct.feature_sigma = feature_sigma;
cluster_struct.clustering_ind = frames;


end

