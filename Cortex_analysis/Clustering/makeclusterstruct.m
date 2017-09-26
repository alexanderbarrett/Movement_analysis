function cluster_struct = makeclusterstruct(cluster_struct_init,labels,feature_mu,feature_sigma,fr,feat_pcs,number_of_clust,feature_labels,clustering_ind,clipped_index_agg)

cluster_struct = cluster_struct_init;
    cluster_struct.labels = labels;
    cluster_struct.feature_mu = feature_mu;
    cluster_struct.feature_sigma = feature_sigma;
    cluster_struct.fr = fr;
    cluster_struct.feat_pcs = feat_pcs;
    cluster_struct.num_clusters = number_of_clust;
    cluster_struct.feature_labels = feature_labels;
        cluster_struct.clustering_ind = clustering_ind;
                cluster_struct.clipped_index_agg = clipped_index_agg;
end