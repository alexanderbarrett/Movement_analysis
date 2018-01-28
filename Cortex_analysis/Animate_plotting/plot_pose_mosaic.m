    
function M = plot_pose_mosaic(mocapstruct,ind_cell)

    %num2str(cluster_numbers');
    fighand =   figure(388);
    set(fighand,'Color','k')
    set(fighand,'Position',[100 100 1100 1100])
    nrows = 8;
    ncols = 8;
    num_vids = min(numel(ind_cell),nrows*ncols);
    
%      for ll = 1:num_vids  %numel(cluster_numbers)
%           %  mov_ind = cluster_numbers(ll);
%             h{ll}=subplot_tight(nrows,ncols, ll);
%      end
%      
    for lk = 1 %numel(movie_output{jjj})
        for ll = 1:num_vids  %numel(cluster_numbers)
          %  mov_ind = cluster_numbers(ll);
           h{ll}=subplot_tight(nrows,ncols, ll);
            movie_size = numel(ind_cell{ll});
            frame_use = mod(lk,movie_size);
            if (frame_use==0)
                frame_use = 1;
            end
            if (numel(ind_cell{ll}))
            plot_frame(mocapstruct,ind_cell{ll}(frame_use),h{ll});
            else
                    set(h{ll},'Color','k')
            end
            %animate_markers_aligned_fullmovie(mocapstruct_all,indshere(1:10:end));

            %imshow(movie_output{mov_ind}(frame_use).cdata,movie_output{mov_ind}(frame_use).colormap)
            if (lk == 1)
            title(num2str(ll),'Color','w')
            end
        end
        M(lk) =      getframe(gcf);
        
    end
    
    %v = VideoWriter(strcat(savedirectory_subcluster,'aggregate_movie_2'),'MPEG-4');
%     open(v)
%     writeVideo(v, M_here)
%     close(v)