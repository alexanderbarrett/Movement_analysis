function examine_features(fighandle,zValues,subset_of_points_to_plot,cond_inds,good_tracks,mocapstructs,videoflag)
  [x, y] = getline(fighandle);
  IN = inpolygon(zValues(:,1),zValues(:,2),...
        x ,y);
    
    beh_list = cell(1,numel(mocapstructs));
    for kk = 1:(numel(mocapstructs))
    [~,ia] = intersect(find(cond_inds==kk),find(IN));
    ia = subset_of_points_to_plot{kk}(ia);
     beh_list{kk} = sort(reshape(unique(rectify_inds(bsxfun(@plus,good_tracks{kk}(ia)',-20:20),max(good_tracks{kk}))),1,[]),'ascend');    
    end
    
    h3=figure(390);
    close 390
    
    
h3=figure(390);
    set(h3,'Color','k')

h1 = subplot(1,2,1)
                    animate_markers_timelapse(mocapstructs{1},beh_list{1}(1:10:end),h1);
                    h2 = subplot(1,2,2)
                    animate_markers_timelapse(mocapstructs{2},beh_list{2}(1:10:end),h2);

    if (videoflag)
               h3=figure(370);
     
          animate_markers_aligned_fullmovie(mocapstructs{1},beh_list{1}(1:10:min(2000,numel(beh_list{1}))));

          animate_markers_aligned_fullmovie(mocapstructs{2},beh_list{2}(1:10:min(2000,numel(beh_list{2}))));
    end
end