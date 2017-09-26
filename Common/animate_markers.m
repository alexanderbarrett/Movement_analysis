function M = animate_markers(markers_preproc,frame_inds,marker_names,markercolor,links,movie_in,save_movie) 

          h=  figure(370)
          set(h,'Color','k')
            plot3( squeeze(markers_preproc.(marker_names{1})(1,1)),...
                squeeze(markers_preproc.(marker_names{1})(1,2)),...
                squeeze(markers_preproc.(marker_names{1})(1,3)),'o','Color',markercolor{1},'MarkerFaceColor',markercolor{1},'MarkerSize',6)
            ax = gca;
            axis(ax,'manual')
            set(gca,'Color','k')
            grid on;
set(gca,'Xcolor',[1 1 1 ]);
set(gca,'Ycolor',[1 1 1]);
set(gca,'Zcolor',[1 1 1]);

            zlim([-40 250])
            xlim([-350 350])
            ylim([-350 350])
            set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
            %matlab_fr = 10;
            
           % frame_inds = time_ordering_fulltrace{jjj}(1:matlab_fr:min(frames_use,numel(time_ordering_fulltrace{jjj})));
            
            %M{jjj} = movie;
            %Mhere = movie;
            frame_last = 0;
            
            marker_plot = ones(1,numel(marker_names));
            
            for lk = frame_inds'%1:10:10000
                
                
                ind_to_plot = lk;
                if (ind_to_plot-frame_last>20)
                    marker_sign = 1;
                end
                
                set(gca,'Nextplot','ReplaceChildren');
                handles_here = cell(1,numel(marker_names));
                for jj = 1:numel(marker_names)
                    % don't plot markers that drop out
                    if (~sum(markers_preproc.(marker_names{jj})(ind_to_plot,:),2) == 0)
                    handles_here{jj} = plot3( squeeze(markers_preproc.(marker_names{jj})(ind_to_plot,1)),...
                        squeeze(markers_preproc.(marker_names{jj})(ind_to_plot,2)),...
                        squeeze(markers_preproc.(marker_names{jj})(ind_to_plot,3)),'o','Color',markercolor{jj},'MarkerFaceColor',markercolor{jj},'MarkerSize',8);
                    hold on
                    marker_plot(jj) = 1;
                    else
                        marker_plot(jj) = 0;
                    end
                end
                
                
                for mm = 1:numel(links)
                    if (ismember(links{mm}(1),1:numel(marker_names)) && ismember(links{mm}(2),1:numel(marker_names)))
                        if (marker_plot(links{mm}(1)) == 1 && marker_plot(links{mm}(2)) == 1)
                    plot3( [squeeze(markers_preproc.(marker_names{links{mm}(1)})(ind_to_plot,1)) ...
                        squeeze(markers_preproc.(marker_names{links{mm}(2)})(ind_to_plot,1)) ],...
                        [ squeeze(markers_preproc.(marker_names{links{mm}(1)})(ind_to_plot,2))...
                        squeeze(markers_preproc.(marker_names{links{mm}(2)})(ind_to_plot,2))],...
                        [squeeze(markers_preproc.(marker_names{links{mm}(1)})(ind_to_plot,3))...
                        squeeze(markers_preproc.(marker_names{links{mm}(2)})(ind_to_plot,3))],'Color',markercolor{links{mm}(1)},'Linewidth',3);
                        end
                    end
                end
                %delete
                if (save_movie)
                    M(find(frame_inds == lk)) = getframe(gcf);
                end
                
                drawnow
                hold off
                
                if (save_movie)
                                    title(strcat('Cluster number: ',num2str(save_movie),'  Frame: ' ,num2str(lk)),'Color','w')
                else
                title(strcat('Frame: ' ,num2str(lk)),'Color','w')
                end
                
                frame_last = lk;
            end
            set(gca,'Nextplot','add');
            
        end