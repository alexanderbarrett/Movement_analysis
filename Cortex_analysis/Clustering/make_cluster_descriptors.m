function make_cluster_descriptors(cluster_struct,feature_plot,savedirectory_subcluster)



 %% compute cluster metrics in feature and cluster space
    figure(33)
    plot(cluster_struct.labels)
   % xlim([1 5000])
    print('-depsc',strcat(savedirectory_subcluster,'labelsovertime.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'labelsovertime.png'))
   
       [idx,C,sumd,D] = kmeans(cluster_struct.feature_mu(:,:),12,'distance','correlation');
    %[idx,C,sumd,D] = kmeans( mean_wave(good_clusters,:),12,'distance','correlation');
    
   [vals,inds] = sort(idx,'ASCEND');
   inds_labels_kmeans = [];
   for lk = 1:numel(inds)
       inds_labels_kmeans = cat(1,inds_labels_kmeans,find(cluster_struct.labels==inds(lk)));
   end
    
    
    cluster_size_dist = cell(1,cluster_struct.num_clusters);
    for ll = 1:cluster_struct.num_clusters
        inst_label = zeros(1,numel(cluster_struct.clustering_ind));
        inst_label(cluster_struct.labels == ll) = 1;
        pixellist = bwconncomp(inst_label);
        cluster_size_dist{ll} = cellfun(@numel,pixellist.PixelIdxList);        
    end
    
    cluster_time_lengths = cellfun(@nanmean,cluster_size_dist);
    
    figure(122)
    hist(cluster_time_lengths,50)
      print('-depsc',strcat(savedirectory_subcluster,'clusterlength_hist.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'clusterlength_hist.png'))
    
    
    figure(123)
    subplot(1,2,1)
    [n,x] = hist(cluster_struct.labels(1:floor(numel(cluster_struct.labels)./2)),1:cluster_struct.num_clusters);
    semilogy(x,n);
    hold on
    [n2,x2] = hist(cluster_struct.labels(floor(numel(cluster_struct.labels)./2):end),1:cluster_struct.num_clusters);
   semilogy(x2,n2,'r');
    
    xlim([1 cluster_struct.num_clusters])
    
    
    figure(55)
    subplot(2,1,1)
    imagesc(cluster_struct.wtAll')
    caxis([0 0.5])
    subplot(2,1,2)
    plot(feature_plot(cluster_struct.clustering_ind))
    
    [~,sorted_ind] = sort(cluster_struct.labels,'ASCEND');
    
    
    figure(34)
    subplot(2,1,1)
    imagesc(cluster_struct.wtAll(sorted_ind,:)')
 caxis([0 3])
    set(gca,'YTick',1:3*numel(cluster_struct.fr):size(cluster_struct.wtAll,1),'YTickLabels',cluster_struct.feature_labels);
    xlabel('Time (frames at 100 fps)')
    
    subplot(2,1,2)
    plot(cluster_struct.labels(sorted_ind))
    xlim([1 numel(cluster_struct.labels)])
    
    print('-depsc',strcat(savedirectory_subcluster,'averagecluster.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'averagecluster.png'))
    
    figure(92)
    plot(0.05.*feature_plot(cluster_struct.clustering_ind),'r')
    hold on
       % plot(markers_preproc.HeadF(clustering_ind,3),'g')

    plot(0.02*cluster_struct.labels-0.003,'b')
    hold off
    
    
    figure(93)
          plot(0.01.*feature_plot(cluster_struct.clustering_ind(sorted_ind)),'r')
    hold on
     %     plot(0.01.*markers_preproc.HeadF(clustering_ind(sorted_ind),3),'g')
    plot(0.02*cluster_struct.labels(sorted_ind)-0.003,'b')
    hold off
    
            print('-depsc',strcat(savedirectory_subcluster,'clustertime_sorted.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'clustertime_sorted.png'))
    
    figure(94)
          plot(0.01.*feature_plot(cluster_struct.clustering_ind( inds_labels_kmeans)),'r')
    hold on
     %     plot(0.01.*markers_preproc.HeadF(clustering_ind(sorted_ind),3),'g')
    plot(0.02*cluster_struct.labels( inds_labels_kmeans)-0.003,'b')
    hold off
      
    print('-depsc',strcat(savedirectory_subcluster,'clustertime_kmeans.eps'))
    print('-dpng',strcat(savedirectory_subcluster,'clustertime_kmeans.png'))