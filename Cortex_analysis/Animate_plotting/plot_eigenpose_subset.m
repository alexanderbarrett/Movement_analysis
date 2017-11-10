function plot_eigenpose_subset(eigenpose,cluster_marker_inds,markercolor,links,h)
    hold on
    set(h,'Color','k')
    plot3( squeeze(eigenpose(1,1)),...
        squeeze(eigenpose(1,2)),...
        squeeze(eigenpose(1,3)),'o','Color',markercolor{cluster_marker_inds(1)},'MarkerFaceColor',markercolor{cluster_marker_inds(1)},'MarkerSize',6)
    ax = gca;
    axis(ax,'manual')
    set(gca,'Color','k')
    grid on;
    set(gca,'Xcolor',[1 1 1 ]);
    set(gca,'Ycolor',[1 1 1]);
    set(gca,'Zcolor',[1 1 1]);
    
    zlim([-125 100])
    xlim([-100 100])
    ylim([-100 100])
    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
    
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
                squeeze(eigenpose(jj,3)),'o','Color',markercolor{cluster_marker_inds(jj)},...
                'MarkerFaceColor',markercolor{cluster_marker_inds(jj)},'MarkerSize',5);
            hold on
          %  marker_plot(jj) = 1;      
    end
    
    
    for mm = 1:numel(links)
        if (ismember(links{mm}(1),cluster_marker_inds) && ismember(links{mm}(2),cluster_marker_inds))
           % if (marker_plot(links{mm}(1)) == 1 && marker_plot(links{mm}(2)) == 1)
           ind_1 = find(cluster_marker_inds==(links{mm}(1)));
                      ind_2 = find(cluster_marker_inds==(links{mm}(2)));

                plot3( [squeeze(eigenpose(ind_1,1)) ...
                    squeeze(eigenpose(ind_2,1)) ],...
                    [ squeeze(eigenpose(ind_1,2))...
                   squeeze(eigenpose(ind_2,2))],...
                    [squeeze(eigenpose(ind_1,3))...
                    squeeze(eigenpose(ind_2,3))],'Color',markercolor{(links{mm}(1))},'Linewidth',3);
           % end
        end
    end