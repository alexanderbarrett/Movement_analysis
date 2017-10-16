function  plot_tsne(agg_features,num_tsne_timepoints,labels)


time_pts = randsample(1:size(agg_features,1),num_tsne_timepoints);

%mappedX = tsne(agg_features(time_pts,:),labels);
mappedX = tsne(agg_features(time_pts,:));

label_color = parula(max(labels(time_pts)));
metacluster_color = parula(10);

%% feature mean
feature_mean = zeros(size(agg_features,2),max(labels));
good_clusters = [];
for zz = 1:max(labels)
    ind_here = find(labels(time_pts) == zz);
        if numel(ind_here)
feature_mean(:,zz) = mean(agg_features(ind_here,:))';
good_clusters = cat(1,good_clusters,zz);
        end   
end

    [idx,C,sumd,D] = kmeans(bsxfun(@rdivide,bsxfun(@minus,feature_mean(:,good_clusters),mean(feature_mean(:,good_clusters),2)),mean(feature_mean(:,good_clusters),2))',min(10,max(labels)),'distance','correlation');    
    [vals,inds] = sort(idx,'ASCEND');

    sorted_labels = (good_clusters(inds));
    sorted_metaclusterind = idx(inds);
    
    
for zz = 1:numel(sorted_labels)
    ind_here = find(labels(time_pts) == sorted_labels(zz));
    if numel(ind_here)
        
figure(44)
plot(mappedX(ind_here,1),mappedX(ind_here,2),'o','Color',label_color(zz,:),'MarkerFaceColor',label_color(zz,:),'MarkerSize',2)
% plot(mappedX(ind_here,1),mappedX(ind_here,2),'o','Color',...
%     metacluster_color(sorted_metaclusterind((zz)),:),...
% 'MarkerFaceColor',metacluster_color(sorted_metaclusterind((zz)),:),'MarkerSize',2)

hold on
    end
end

end