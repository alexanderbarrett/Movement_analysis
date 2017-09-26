function feature_euclidean_distance(cluster_struct,cluster_struct2,wtAll,wtAll2,savedirectory_subcluster)

  feature_mu_reshape = reshape(cluster_struct.feature_mu,size(cluster_struct.feature_mu,1),numel(cluster_struct.fr),[]);
   feature_mu_reshape2 = reshape(cluster_struct2.feature_mu,size(cluster_struct2.feature_mu,1),numel(cluster_struct2.fr),[]);
   
   agg_feature_mu_reshaped = cat(1,feature_mu_reshape,feature_mu_reshape2);
   agg_feature_mu = cat(1,cluster_struct.feature_mu,cluster_struct2.feature_mu);
   total_clusters = size(agg_feature_mu,1);
  % wtAll_total = cat(1,wtAll,wtAll2);
   
    x_ind = size(feature_mu_reshape,2);
    y_ind = size(feature_mu_reshape,3);
    
    xcorr_features = zeros(total_clusters,total_clusters);
    euc_features = zeros(total_clusters,total_clusters);
    wave_features = zeros(total_clusters,total_clusters);
    mean_wave = zeros(total_clusters,size(wtAll,2));
    for lkl = 1:size(cluster_struct.feature_mu,2)
        mean_wave(lkl,:) = nanmean(wtAll((cluster_struct.labels==lkl),:),1);
    end
    
     for lkl = size(cluster_struct.feature_mu,1)+1:(size(cluster_struct.feature_mu,1)+size(cluster_struct2.feature_mu,1))
        mean_wave(lkl,:) = nanmean(wtAll2(cluster_struct2.labels==(lkl-size(cluster_struct.feature_mu,1)),:),1);
    end
    
    
    for jjj = 1:total_clusters;
        for mkm = jjj:total_clusters
            corr_c = xcorr2(squeeze(agg_feature_mu_reshaped(jjj,:,:)),squeeze(agg_feature_mu_reshaped(mkm,:,:)));
            xcorr_features(jjj,mkm) = corr_c(x_ind,y_ind);
            euc_features(jjj,mkm) = pdist2(agg_feature_mu(jjj,:),agg_feature_mu(mkm,:),'correlation');
            
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

    [idx,C,sumd,D] = kmeans(agg_feature_mu,24,'distance','correlation');
    % [idx,C,sumd,D] = kmeans( mean_wave(good_clusters,:),12,'distance','correlation');
    
    [vals,inds] = sort(idx,'ASCEND');
    
    
    figure(371)
    imagesc(xcorr_features(:,:))
    %caxis([-10 10])
    print('-depsc',strcat(savedirectory_subcluster,'xcorrclustered','.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'xcorrclustered','.png'))
    
    figure(373)
    imagesc(euc_features(:,:))
    caxis([0 3])
    print('-depsc',strcat(savedirectory_subcluster,'eucclustered','.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'eucclustered','.png'))
    
       figure(374)
    imagesc(euc_features(inds,inds))
    caxis([0 3])
    print('-depsc',strcat(savedirectory_subcluster,'eucclustered_sorted','.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'eucclustered_sorted','.png'))
    
    figure(372)
    imagesc(wave_features(:,:))
    caxis([0 0.4])
end