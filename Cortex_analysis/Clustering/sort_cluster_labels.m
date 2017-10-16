function sorted_clusters = sort_cluster_labels(cluster_struct)


good_clusters = cluster_struct.labels;
    
   
    feature_mu_reshape = reshape(cluster_struct.feature_mu,size(cluster_struct.feature_mu,1),numel(cluster_struct.fr),[]);
    
    x_ind = size(feature_mu_reshape,2);
    y_ind = size(feature_mu_reshape,3);
    
    xcorr_features = zeros(numel(good_clusters),numel(good_clusters));
    euc_features = zeros(numel(good_clusters),numel(good_clusters));
    wave_features = zeros(numel(good_clusters),numel(good_clusters));
    mean_wave = zeros(max(good_clusters),size(cluster_struct.wtAll,2));
    for lkl = 1:max(good_clusters)
        mean_wave(lkl,:) = mean(cluster_struct.wtAll(find(cluster_struct.labels==lkl),:),1);
    end
    
    
    for jjj = good_clusters(1:end);
        for mkm = good_clusters(find(good_clusters == jjj):end)
            corr_c = xcorr2(squeeze(feature_mu_reshape(jjj,:,:)),squeeze(feature_mu_reshape(mkm,:,:)));
            xcorr_features(jjj,mkm) = corr_c(x_ind,y_ind);
            euc_features(jjj,mkm) = pdist2(cluster_struct.feature_mu(jjj,:),cluster_struct.feature_mu(mkm,:),'correlation');
            
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

    [idx,C,sumd,D] = kmeans(cluster_struct.feature_mu(good_clusters,:),min(12,numel(good_clusters)),'distance','correlation');
    % [idx,C,sumd,D] = kmeans( mean_wave(good_clusters,:),12,'distance','correlation');
    
    [vals,inds] = sort(idx,'ASCEND');
    
    
    figure(371)
    imagesc(xcorr_features(good_clusters(inds),good_clusters(inds)))
    caxis([-10 10])
    print('-depsc',strcat(savedirectory_subcluster,'xcorrclustered','.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'xcorrclustered','.png'))
    
    figure(373)
    imagesc(euc_features(good_clusters(inds),good_clusters(inds)))
    caxis([0 3])
    print('-depsc',strcat(savedirectory_subcluster,'eucclustered','.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'eucclustered','.png'))
    
    figure(372)
    imagesc(wave_features(good_clusters(inds),good_clusters(inds)))
    caxis([0 0.4])
    
    
sorted_clusters =good_clusters(inds);