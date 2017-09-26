

function [cluster_struct] = Cluster_GMM(agg_features,opts,frames)
%simple GMM clustering
% inputs: 


 %% could do a 'pre-pca' stage. To ignore this, set num_pcs1 high
        [coeff,score,latent,tsquared,explained] = pca(agg_features(:,frames)');
     %   num_pcs_1 = size(agg_features,1);
     
        
        pc_traces = score';
        opts.num_pcs_1 = min(opts.num_pcs_1,size(agg_features,1));
        
        pc_spectrograms = cell(1,opts.num_pcs_1);
        
        if ~size(score,2)
            fprintf('UH OH ERROR IN PCA \n')
        else
        
        for k=1:opts.num_pcs_1
          %  [~,~,time_clustering,pc_spectrograms{k}] = spectrogram(agg_features(k,:),clustering_window,clustering_overlap,fps,fps);
            [~,~,time_clustering,pc_spectrograms{k}] = spectrogram(pc_traces(k,:),opts.clustering_window,opts.clustering_overlap,opts.fps,opts.fps);
            
%         
%             parameters.numPeriods = 50;
%             parameters.samplingFreq = 245;
%             parameters.maxF = 100;
%             parameters.minF = 1;
%             parameters.omega0 = 5;
%             
%             [pc_spectrograms{k},f,scales] = findWavelets(agg_features(k,:)',20,parameters);
%             pc_spectrograms{k} = pc_spectrograms{k}';
        end
        fr = 0:opts.fps./2;
        num_fr = numel(fr);
        agg_spectrograms = cell2mat(pc_spectrograms'); %second dimension is time base
        
        timebase = size(pc_spectrograms{1},2);
        
        frames_per_bin = floor(size(agg_features,2)./timebase);
        
        [coeff2,score2,latent2,tsquared2,explained2] = pca(agg_spectrograms(:,:)');
        pca_agg_spectrogram = score2';
        
        %% do a GMM clustering
        obj = fitgmdist((pca_agg_spectrogram(1:opts.num_pcs_2 ,:))',opts.num_clusters,...
            'Start', 'plus', 'Options', statset('MaxIter',1000),'RegularizationValue',0.1);
        clusterobj = cluster(obj,(pca_agg_spectrogram(1:opts.num_pcs_2 ,:))');
        
                labels =clusterobj;
                        featurespectrogram = pca_agg_spectrogram(1:opts.num_pcs_2 ,:);


                        feat_pcs_full = reshape(coeff2,num_fr,opts.num_pcs_1,size(coeff2,2));
     cluster_struct.feat_pcs = feat_pcs_full(:,:,1:opts.num_pcs_2);
% feat_use = permute(cluster_struct.feat_pcs,[3 1 2]);
%                         
%                         %% get mu and var
%                         mu = obj.mu;
% sigma = obj.Sigma;
% feature_mu = mtimesx(mu,feat_use);
% 
% figure(44)
% imagesc(log10(abs(squeeze(feature_mu(6,:,:))))')
% 
% feature_sigma = zeros(size(feat_use,2),size(feat_use,3) , numclusters,size(feat_use,2),size(feat_use,3));
% for j = 1:numclusters
% feature_sigma(:,:,j,:,:) = mtimesx(permute(feat_use,[2 3 1]),mtimesx(squeeze(sigma(:,:,j),feat_use)));
% end
%                         
%         
V = coeff2(1:opts.num_pcs_2);

mu = obj.mu;
sigma = obj.Sigma;

feature_mu = mu * V';
feature_sigma = zeros(size(featurespectrogram,1),size(featurespectrogram,1) ,opts.num_clusters);
for j = 1:opts.num_clusters
feature_sigma(:,:,j) = V * sigma(:,:,j) * V';
end



                        
cluster_struct = struct();
    cluster_struct.labels = labels;
        cluster_struct.num_clusters = opts.num_pcs_2;
cluster_struct.wtAll = featurespectrogram';
     cluster_struct.fr = fr;
     
     cluster_struct.feature_mu = feature_mu;
         cluster_struct.feature_mu = feature_mu;
    cluster_struct.feature_sigma = feature_sigma;
     cluster_struct.clustering_ind = frames;

     
%     cluster_struct.feature_labels = feature_labels;
%         cluster_struct.clustering_ind = 
%                 cluster_struct.clipped_index_agg = clipped_index_agg;
                        
        end
end
        