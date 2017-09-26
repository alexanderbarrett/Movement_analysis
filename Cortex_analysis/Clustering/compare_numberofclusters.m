function hyperparam_comparison_struct= compare_numberofclusters(cluster_struct,modular_cluster_properties,modno,cluster_scan,opts)%agg_features,cluster_scan)

    num_iter = numel(cluster_scan);
    shortCount_agg = zeros(1,num_iter);
  %  opts.num = size(modular_cluster_properties.agg_features{modno},1); % number modes (spectrograms to find) (usually want full dimension)
    
    
    markovLLR_agg = zeros(1,num_iter);
    entropy_agg = zeros(1,num_iter);
    exits_agg = zeros(1,num_iter);
    meanDwell_agg = zeros(1,num_iter);
    meanClusterDist_agg = zeros(1,num_iter);
    cluster_time_lengths_agg = cell(1,num_iter);
    fano_agg =cell(1,num_iter);
    feature_similarity_agg =cell(1,num_iter);
    feature_similarity_summary = cell(1,num_iter);
    feature_similarity_summary_tot = zeros(1,num_iter);
    
        
        for nk = 1:num_iter
            opts.numclusters =   cluster_scan(nk);
                        fr = cluster_struct.fr;
                        wtAll = cluster_struct.wtAll;
clustering_ind = cluster_struct.clustering_ind;
%             %% make
%             if (nk == 1)
%                 [labels , feature_mu, feature_sigma,wtAll,fr]= WaveletCluster(agg_features(:,clustering_ind)',1:numel(clustering_ind),opts);
%             else
                [cluster_struct_int]= WaveletCluster(modular_cluster_properties.agg_features{modno}(:,clustering_ind)',1:numel(clustering_ind),opts,wtAll,fr);
                   
                labels = cluster_struct_int.labels;
%                 
%                 cluster_struct.labels = labels;
%     cluster_struct.feature_mu = feature_mu;
%     cluster_struct.feature_sigma = feature_sigma;
%     cluster_struct.fr = fr;
%     cluster_struct.wtAll = wtAll;
           % end
            
            
            
            cluster_size_dist = cell(1,cluster_scan(nk));
            for ll = 1:cluster_scan(nk)
                inst_label = zeros(1,numel(clustering_ind));
                inst_label(labels == ll) = 1;
                pixellist = bwconncomp(inst_label);
                cluster_size_dist{ll} = cellfun(@numel,pixellist.PixelIdxList);
            end
            
            cluster_time_lengths_agg{nk} = cellfun(@nanmean,cluster_size_dist);
            emptyclusters = find(isnan(cluster_time_lengths_agg{nk})==1);
            label_iter = 0;
            for mm = emptyclusters
                labels(find(labels>(mm-label_iter))) = labels(find(labels>(mm-label_iter)))-1;
                label_iter = label_iter+1;
            end
            
            metrics_struct = clusterMetricsTodd_JDM(cluster_struct_int);
            
            shortCount_agg(nk) =  metrics_struct.shortCount;
                markovLLR_agg(nk)=  metrics_struct.markovLLR;
                entropy_agg(nk)=  metrics_struct.entropy;
                exits_agg(nk)=  metrics_struct.exits;
                meanDwell_agg(nk)=  metrics_struct.meanDwell;
              meanClusterDist_agg(nk)=  metrics_struct.meanClusterDist;
            
            
            
            feature_fano = zeros(1,cluster_scan(nk));
            feature_similarity_score = zeros(cluster_scan(nk),cluster_scan(nk)); %(jk,ll)
            for jk = 1:cluster_scan(nk)
                cluster_mean = mean(wtAll(find(labels == jk),:),1);
                cluster_std = std(wtAll(find(labels == jk),:),[],1);
                
                feature_fano(jk) = mean(cluster_std./abs(cluster_mean));
                for ll=jk:cluster_scan(nk)
                    second_index = find(labels==ll);
                    feature_similarity_score(jk,ll) = mean(mean(pdist2(wtAll(second_index(1:10:end),:),cluster_mean,...
                        'euclidean')));
                end
            end
            
            feature_similarity_agg{nk} = feature_similarity_score;
            
            fano_agg{nk} = feature_fano;
            end
            
            for nk = 1:num_iter
                
                for jk = 1:cluster_scan(nk)
                    offdiag = cat(2,feature_similarity_agg{nk}(jk,:),feature_similarity_agg{nk} (:,jk)');
                    offdiag(offdiag == 0) = [];
                    feature_similarity_summary{nk}(jk,1) = feature_similarity_agg{nk} (jk,jk);
                    feature_similarity_summary{nk}(jk,2) = nanmean(offdiag);
                end
                feature_similarity_summary_tot(nk) =...
                    nanmean(feature_similarity_summary{nk}(:,1)./feature_similarity_summary{nk}(:,2));
            end
        

        
        hyperparam_comparison_struct.markovLLR_agg =markovLLR_agg;
        hyperparam_comparison_struct.entropy_agg = entropy_agg;
        hyperparam_comparison_struct.exits_agg = exits_agg;
        hyperparam_comparison_struct.meanDwell_agg = meanDwell_agg;
        hyperparam_comparison_struct.meanClusterDist_agg = meanClusterDist_agg;
        hyperparam_comparison_struct.cluster_time_lengths_agg = cluster_time_lengths_agg;
       hyperparam_comparison_struct.fano_agg = fano_agg;
        hyperparam_comparison_struct.feature_similarity_agg = feature_similarity_agg;
        hyperparam_comparison_struct.feature_similarity_summary = feature_similarity_summary;
           hyperparam_comparison_struct.feature_similarity_summary_tot = feature_similarity_summary_tot;
end