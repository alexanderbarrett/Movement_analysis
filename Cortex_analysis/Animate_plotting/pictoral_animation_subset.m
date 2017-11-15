function pictoral_animation_subset(agg_features,mocapstruct,frames_ind,offset,cluster_marker_inds,h)


eigenpose = squeeze(agg_features(1,:,:));
subplot(2,1,1)
 hold on
    set(h,'Color','k')
    plot3( squeeze(eigenpose(1,1)),...
        squeeze(eigenpose(1,2)),...
        squeeze(eigenpose(1,3)),'o','Color',mocapstruct.markercolor{cluster_marker_inds(1)},'MarkerFaceColor',mocapstruct.markercolor{cluster_marker_inds(1)},'MarkerSize',6)
    ax = gca;
    axis(ax,'manual')
    set(gca,'Color','k')
    grid on;
    set(gca,'Xcolor',[1 1 1 ]);
    set(gca,'Ycolor',[1 1 1]);
    set(gca,'Zcolor',[1 1 1]);
    
    zlim([-125 100])
    xlim([-100 400])
    ylim([-100 100])
    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])

    
for mm = frames_ind
eigenpose = squeeze(agg_features(mm,:,:));

eigenpose(:,1) = eigenpose(:,1) + offset*(find(frames_ind == mm));

%     
%     ind_to_plot = lk;
%     if (ind_to_plot-frame_last>20)
%         marker_sign = 1;
%     end
    
    %set(gca,'Nextplot','ReplaceChildren');
    handles_here = cell(1,numel(cluster_marker_inds));
    for jj = 1:numel(cluster_marker_inds)
        % don't plot markers that drop out
            handles_here{jj} = plot3( squeeze(eigenpose(jj,1)),...
                squeeze(eigenpose(jj,2)),...
                squeeze(eigenpose(jj,3)),'o','Color',mocapstruct.markercolor{cluster_marker_inds(jj)},...
                'MarkerFaceColor',mocapstruct.markercolor{cluster_marker_inds(jj)},'MarkerSize',5);
            hold on
          %  marker_plot(jj) = 1;      
    end
    
    
    for mm = 1:numel(mocapstruct.links)
        if (ismember(mocapstruct.links{mm}(1),cluster_marker_inds) && ismember(mocapstruct.links{mm}(2),cluster_marker_inds))
           % if (marker_plot(links{mm}(1)) == 1 && marker_plot(links{mm}(2)) == 1)
           ind_1 = find(cluster_marker_inds==(mocapstruct.links{mm}(1)));
                      ind_2 = find(cluster_marker_inds==(mocapstruct.links{mm}(2)));

                plot3( [squeeze(eigenpose(ind_1,1)) ...
                    squeeze(eigenpose(ind_2,1)) ],...
                    [ squeeze(eigenpose(ind_1,2))...
                   squeeze(eigenpose(ind_2,2))],...
                    [squeeze(eigenpose(ind_1,3))...
                    squeeze(eigenpose(ind_2,3))],'Color',mocapstruct.markercolor{(mocapstruct.links{mm}(1))},'Linewidth',3);
           % end
        end
    end
    
end


subplot(2,1,2)
plot(agg_features(frames_ind(1):frames_ind(end),:,3))

end