%% function to coembed two datasets
%s
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';

%% PARMAETERS
%% set parameters
%V8
ratname = 'JDM25';
conditionnumbers = [1,2]; %pre/post
%conditionnumbers = [4,12]; %pre/post task

ratname = 'JDM33'
conditionnumbers = [1,5,10]; %pre/post

ratname = 'Vicon8';
conditionnames = {'Vicon8_caff','Vicon8_prelesion','Vicon8_amph'};
conditionnumbers = [6 1 7];

conditions_to_run = [1:2];
mocapstruct_array = cell(1,numel(conditions_to_run));
ML_features = cell(1,numel(conditions_to_run));

%for referencing after aggregation

move_flag = 1;


overwrite_MLfile = 1;
overwrite_MLcoeff_run = 0;
%% Load files
for ll =1:numel(conditions_to_run)
    % overwrite coeffs on the first file
    if overwrite_MLcoeff_run
    overwrite_MLfile = 1;
    if ll == 1
        overwrite_MLcoeff = 1;
    else
        overwrite_MLcoeff = 0;
    end
    else
    overwrite_MLfile = 0;
            overwrite_MLcoeff = 0;
    end
    
    descriptor_struct = get_mocap_files_table(conditionnumbers(conditions_to_run(ll)),ratname);
    [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
    [mocapstruct_array{ll}] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,0,[],mocapvideodirectory,1);
    %% load in the ML features
    ML_features{ll} = get_supervised_features(mocapstruct_array{ll},...
        mocapstruct_array{ll}.modular_cluster_properties.clustering_inds_agg{2},2,ratname,'temp',overwrite_MLfile,overwrite_MLcoeff);
end
    
tsnegranularity = 15;

%recompute these 
condition_inds =  [];
frames_with_good_tracking = cell(1,numel(conditions_to_run));
ml_subset_moving = cell(1,numel(conditions_to_run));
candidate_features = cell(1,numel(conditions_to_run));
jt_features = [];

%% set tsne parameters
subset_of_points_to_plot_tsne = 1:tsnegranularity:100*25000;

[ tsne_features,num_feat]= load_tsne_features('annotation_ii');
numDims = 2; pcaDims = min(30,sum(num_feat)+numel(tsne_features)-numel(num_feat));
perplexity = 30; theta = .5; alg = 'svd';

for ll =1:numel(conditions_to_run)
    tsnefeat_name = cell(0,1);

    frames_with_good_tracking{ll} = ML_features{ll}.framelist_true;
    %% aggregate features
    if ~move_flag
        for mm = 1:numel(tsne_features)
              if mm>numel(num_feat)
          num_feat(mm) = 1;
      end
            candidate_features{ll} = cat(2,candidate_features{ll} ,candidate_features{ll}.(tsne_features{mm})(subset_of_points_to_plot_tsne,1:num_feat(mm)));
            ml_subset_moving{ll} = 1:size(candidate_features{ll},1);
        end
    else
        [frames_with_good_tracking{ll},~,ml_subset_moving{ll} ] = intersect(mocapstruct_array{ll}.move_frames,frames_with_good_tracking{ll});
        for mm = 1:numel(tsne_features)
              if mm>numel(num_feat)
          num_feat(mm) = 1;
              end
      
            candidate_features{ll} = cat(2,candidate_features{ll} ,ML_features{ll}.(tsne_features{mm})(ml_subset_moving{ll}(subset_of_points_to_plot_tsne),1:num_feat(mm)));
            for jj = 1:num_feat(mm)
            tsnefeat_name{end+1} = strcat(tsne_features{mm},'__',num2str(jj));
            end
        end
    end
    
    %whiten features
    candidate_features{ll} =bsxfun(@minus,candidate_features{ll},...
      nanmean(candidate_features{ll},1));
    %aggregate features
    jt_features = cat(1,jt_features,(candidate_features{ll}));
    
    %good_track_switch = {frames_with_goodtracking_caff,[],frames_with_goodtracking_amph };
condition_inds = cat(1,condition_inds,ll*ones(numel(subset_of_points_to_plot_tsne),1));
end

subset_of_points_to_plot_tsne_move = {ml_subset_moving{1}(subset_of_points_to_plot_tsne)};% ml_subset_moving{2}(subset_of_points_to_plot_tsne)};

% whiten with the std here 
  jt_features = bsxfun(@rdivide,bsxfun(@minus,jt_features,...
      nanmean(jt_features,1)),nanstd(jt_features,[],1));

  figure(98)
  plot(jt_features(:,2))
  
%% start Tsne
fprintf('starting tsne for %f frames  \n',size(jt_features,1));
tic
zValues =  fast_tsne(jt_features, numDims, pcaDims, perplexity, theta, alg);
toc

zvalues_1_only = fast_tsne(bsxfun(@rdivide,candidate_features{1},nanstd(candidate_features{1})), numDims, pcaDims, perplexity, theta, alg);

h2=figure(608)
plot(zvalues_1_only(:,1),zvalues_1_only(:,2),'b+')

%% plot the results
h=figure(605)
plot(zValues(find(condition_inds==1),1),zValues(find(condition_inds==1),2),'b+')
hold on
plot(zValues(find(condition_inds==2),1),...
    zValues(find(condition_inds==2),2),'r+')
hold off


beh_list = examine_features(h,zValues,subset_of_points_to_plot_tsne_move,...
    condition_inds,frames_with_good_tracking,mocapstruct_array,1);


   

examine_clusters = 0;
if examine_clusters
[x, y] = getline(h);
IN = inpolygon(zValues(:,1),zValues(:,2),...
    x ,y);


beh_list = cell(1,3);
for kk = [1 3]
    [~,ia] = intersect(find(condition_inds==kk),find(IN));
    ia = subset_of_points_to_plot_tsne(ia);
    beh_list{kk} = sort(reshape(unique(rectify_inds(bsxfun(@plus,good_track_switch{kk}(ia)',-10:10),max(good_track_switch{kk}))),1,[]),'ascend');
end
% close 370
h3=figure(370)
set(h3,'Color','k')

h1 = subplot(1,2,1)
animate_markers_timelapse(mocapstruct_caff,beh_list{1}(1:10:end),h1);
h2 = subplot(1,2,2)
animate_markers_timelapse(mocapstruct_prel,beh_list{2}(1:10:end),h2);
animate_markers_timelapse(mocapstruct_amph,beh_list{3}(1:10:end),h2);


animate_markers_aligned_fullmovie(mocapstruct_prel,beh_list{2}(1:10:end));
animate_markers_aligned_fullmovie(mocapstruct_caff,beh_list{1}(1:10:end));
animate_markers_aligned_fullmovie(mocapstruct_amph,beh_list{3}(1:10:end));
end


%% look at embedded space velocity (optional()
%examine_embedded_velocity(zValues)


%% density, watershed, images

tsnehere =zvalues_1_only;
%tsnehere = zValues;

density_maps = cell(1,numel(conditions_to_run));
density_res = 1001;
density_max = max(tsnehere(:));
density_width = 0.75;

conditions_to_run = 1;
figure(480)
for ll =1:numel(conditions_to_run)
[xx,yy,density_maps{ll}] = findPointDensity(tsnehere(find(condition_inds==ll),:),density_width,[density_res density_res],[-density_max density_max]);
subplot(1,numel(conditions_to_run)+1,ll)
imagesc(flipud(density_maps{ll}))
%to try [f,xi] = ksdensity(x)
end
[xx,yy,density_jt] = findPointDensity(tsnehere(:,:),density_width,[density_res density_res],[-density_max density_max]);
subplot(1,numel(conditions_to_run)+1,numel(conditions_to_run)+1)
imagesc(flipud(density_jt))




%% now threshold the density maps and apply the watershed
figure(482)
[x,n] = hist(reshape((density_maps{2}),1,[]),1000);
loglog(n,x)
figure(582)
[x,n] = hist(reshape((density_maps{1}),1,[]),logspace(-5,-3,200));
[x2,n2] = hist(reshape((density_maps{2}),1,[]),logspace(-5,-3,200));

loglog(n,x)
hold on
loglog(n2,x2,'r')
hold off

density_thresholds{1} = 3*10^(-5);
density_thresholds{2} = 3*10^(-5);


%figure(490)
num_conditions = numel(conditions_to_run);
pdf_maps = cell(1,num_conditions);
density_watersheds = cell(1,num_conditions);
density_cc = cell(1,num_conditions);
density_stats = cell(1,num_conditions);
for ll = [1 num_conditions]
    density_maps{ll}(density_maps{ll}<density_thresholds{ll}) = 0;
    inv_density_jt = 1./density_maps{ll};
    inv_density_jt(inv_density_jt>(1/density_thresholds{ll})) = 10^10;
    L=watershed(inv_density_jt);
    L(density_maps{ll}==0) = 0;
    density_cc{ll} = bwconncomp(L);
    density_stats{ll} = regionprops(density_cc{ll},'convexhull');
    density_watersheds{ll} = L;
end

figure(482)
subplot(1,3,1)
imagesc(flipud(density_watersheds{1}))
subplot(1,3,2)
imagesc(flipud(density_watersheds{2}))

L = flipud(density_watersheds{cond_select});
s = regionprops(L, 'Centroid');

for k = 1:numel(s)
    c = s(k).Centroid;
    text(c(1), c(2), num2str(k), ... %sprintf('%d', integer(ind)),
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','white');
end
hold off

%% clustering
cond_select = 1;
density_objects = density_cc{cond_select}.NumObjects;

example_inds = cell(1,density_objects);
cluster_num = zeros(1,density_objects);
cluster_mean_length = zeros(4,density_objects);
example_inds_cond = cell(density_objects,num_conditions);
number_cond = zeros(density_objects,num_conditions+1);
example_inds_cond_video = cell(density_objects,num_conditions);
annotation_vec = cell(1,num_conditions);
annotation_labels = cell(1,density_objects);
annotation_labels{1} = 'null';

for ll = 1:num_conditions
    annotation_vec{ll} = zeros(1,max(subset_of_points_to_plot_tsne_move{ll}));
end


subset_of_points_to_plot =1:size(tsnehere,1);
condition_labels = cat(1,ones(1,numel(subset_of_points_to_plot)));%,2*ones(1,numel(subset_of_points_to_plot)));

for jj = 1:density_objects
    [example_inds_x,example_inds_y]  = ind2sub( density_cc{cond_select}.ImageSize,...
        density_cc{cond_select}.PixelIdxList{jj} );
    verts_x = (round(density_stats{cond_select}(jj).ConvexHull(:,1)));
    verts_x(verts_x<1) = 1;
    verts_x(verts_x>numel(xx)) =numel(xx);
    
    verts_y = (round(density_stats{cond_select}(jj).ConvexHull(:,2)));
    verts_y(verts_y<1) = 1;
    verts_y(verts_y>numel(yy)) =numel(yy);
    %% get the points in the watershed polygon
    IN = inpolygon(tsnehere(:,1),tsnehere(:,2),...
        xx(verts_x) ,yy(verts_y));
    
    % IN = intersect(find(IN), find(density_jt >0));
    %% get the bout characteristics
    example_inds{jj} =  find(IN);
    annotation_labels{jj+1} = strcat('cluster ',num2str(jj));
    
    for ll = 1:num_conditions
        [~,example_inds_cond{jj,ll}] = intersect(find(condition_inds==ll),example_inds{jj});
        annotation_vec{ll}(subset_of_points_to_plot_tsne_move{ll}(example_inds_cond{jj,ll}))=jj;
    end
    %number_cond(jj,:) = [numel(example_inds_cond{jj,1}) numel(example_inds_cond{jj,2})  jj];
    
    
    
    
    if numel(example_inds{jj})
        temp = rectify_inds(unique(bsxfun(@plus,example_inds{jj},-10:10)),max(subset_of_points_to_plot_tsne));
        test = zeros(1,max(subset_of_points_to_plot_tsne));
        test(temp) = 1;
        characteristics = bwconncomp(test);
        cluster_num(jj) = characteristics.NumObjects;
        cluster_mean_length(:,jj) = [characteristics.NumObjects mean(cellfun(@numel,characteristics.PixelIdxList)) ...
            median(cellfun(@numel,characteristics.PixelIdxList)) std(cellfun(@numel,characteristics.PixelIdxList))];
        
        for kk =1:num_conditions
            if numel(example_inds_cond{jj,kk})
            thresh = 2;
            diff_input = cat(1,0,diff(example_inds_cond{jj,kk}));
            diff_input(diff_input>thresh) = 0;
            outputstruct = bwconncomp(diff_input);
           % indivbouts = cell(1,numel(outputstruct.PixelIdxList));
            fulloutput = [];
            for mm = 1:numel(outputstruct.PixelIdxList)
                fulloutput = cat(2,fulloutput,...
                    frames_with_good_tracking{kk}(subset_of_points_to_plot_tsne_move{kk}(example_inds_cond{jj,kk}(outputstruct.PixelIdxList{mm}(1)))):...
                    frames_with_good_tracking{kk}(subset_of_points_to_plot_tsne_move{kk}(example_inds_cond{jj,kk}(outputstruct.PixelIdxList{mm}(end)))));                
            end
            %pick up single instances
              fulloutput = cat(2,fulloutput,frames_with_good_tracking{kk}(subset_of_points_to_plot_tsne_move{kk}(example_inds_cond{jj,kk}(diff_input == 0))));
            
            temp = (unique(bsxfun(@plus,fulloutput',-20:20)));
            temp(temp<1) = 1;
            example_inds_cond_video{jj,kk} = temp;
            else
                example_inds_cond_video{jj,kk} =[];
            end
        end
    end
end
good_clusters = find(cellfun(@numel,example_inds_cond(:,cond_select))>100)';
%sum()
make_dotplot(annotation_vec{cond_select},annotation_labels,ML_features{cond_select},good_clusters)
annotation_vec{cond_select}(annotation_vec{cond_select}==0) = [];
annot_reordered = reorder_annotation_vec(annotation_vec{cond_select},sorted_clust_ind);
compute_transition_matrix(annot_reordered,annotation_labels)

animate_markers_aligned_fullmovie(mocapstruct_caff,example_inds_cond_video{170,1}(1:10:end));
animate_markers_aligned_fullmovie(mocapstruct_amph,example_inds_cond_video{170,3}(1:10:end));


%compute_behavior_variances
var_array = zeros(10,density_objects);
mean_pose = zeros(36,density_objects);
pose_var = zeros(36,density_objects);
frames_agg =cell(1,density_objects);
for ll = 1:density_objects %good_clusters %
    mlist = [1:10,17,18];
    allframes = [];
    for jj = mlist
        allframes =  cat(2,allframes,...
            mocapstruct_array{cond_select}.markers_aligned_preproc.(mocapstruct_array{cond_select}.markernames{jj})(example_inds_cond_video{ll,cond_select},:));
    end
    %   frames_agg{ll} = allframes;
    mean_pose(:,ll) = mean(allframes,1);
    pose_var(:,ll) = std(allframes,[],1);
    
    %   [COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(allframes);
    %  var_array(:,ll) =  (EXPLAINED(1:10))';
end
pose_dist = squareform(pdist(mean_pose','euclidean'));
figure(1112)
Z = linkage(mean_pose','ward','euclidean');
c = cluster(Z,'maxclust',15);
[newind,sorted_clust_ind] = sort(c,'ASCEND');
%crosstab(c,species)
[H,T,outperm] = dendrogram(Z,density_objects);
xlabel('cluster number (original)')

figure(1111)
imagesc(pose_dist(outperm,outperm))
imagesc(pose_dist(sorted_clust_ind,sorted_clust_ind))

caxis([0 150])
sorted_clust_ind= outperm;

L = zeros(size(density_watersheds{cond_select}));
Lnew = zeros(size(L));
Lnew2 = (zeros(size(L,1),size(L,2),3));
hold on
vv = parula(93);
for k = 1:numel(sorted_clust_ind)
    L(density_cc{cond_select}.PixelIdxList{k}) = k;
end
for k = 1:numel(sorted_clust_ind)
    indstorelabel = find(L ==  sorted_clust_ind(k));%density_cc{cond_select}.PixelIdxList{sorted_clust_ind(k)};
    Lnew(indstorelabel) = (k);
   % [a,b] = find(Lnew == k);
    %% all i can say is fuck matlab
    for zz = 1:numel(a)
   % Lnew2(a(zz),b(zz),:) =(vv(k,:));
    end
end

figure(487)
subplot(1,2,1)
im=imagesc(flipud(Lnew))
colorbar

subplot(1,2,2)
imagesc(flipud(L))
colorbar

ind = 0;
s = regionprops(flipud(L), 'Centroid');
s2 = regionprops(flipud(Lnew), 'Centroid');

for k = 1:numel(sorted_clust_ind')
    
    % get the ind of the sorted cluster
            c = s2((k)).Centroid;

    subplot(1,2,1)
    text(c(1), c(2), num2str((k)), ... %sprintf('%d', integer(ind)),
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','white');
    
        c = s((k)).Centroid;

        % give it the sorted label
    subplot(1,2,2)
    text(c(1), c(2), num2str((k)), ... %sprintf('%d', integer(ind)),
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','white');
end
hold off
textColor = 'white';


figure(89)
hold on

for jj = 1:density_objects
    [example_inds_x,example_inds_y]  = ind2sub( density_cc{cond_select}.ImageSize,...
        density_cc{cond_select}.PixelIdxList{jj} );
    verts_x = (round(density_stats{cond_select}(jj).ConvexHull(:,1)));
    verts_x(verts_x<1) = 1;
  verts_x(verts_x>numel(xx)) =numel(xx);
    
    verts_y = (round(density_stats{cond_select}(jj).ConvexHull(:,2)));
    verts_y(verts_y<1) = 1;
    verts_y(verts_y>numel(yy)) =numel(yy);
    plot(xx(verts_x),yy(verts_y),'k')
end
plot(tsnehere(find(condition_inds==1),1),tsnehere(find(condition_inds==1),2),'b+')

hold off

hold on
%plot(tsnehere(find(condition_inds==2),1),...
%    tsnehere(find(condition_inds==2),2),'r+')
%s2 = regionprops(Lnew, 'Centroid');
s2 = regionprops(Lnew, 'Centroid');

for k = 1:numel(sorted_clust_ind')
    
    % get the ind of the sorted cluster
            c = s2((k)).Centroid;

    text(xx(floor(c(1))), yy(floor(c(2))), num2str((k)), ... %sprintf('%d', integer(ind)),
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','black');
end
hold off



showvideo = 1;
plot_timelapse_mosaic(mocapstruct_array{cond_select},example_inds_cond_video(sorted_clust_ind(1:end),1))
videodir = 'Y:\Jesse\Data\Motionanalysis_captures\test_clustervideos';
videotag = 'caff_solo_embed_dense';
save_clustered_movies(mocapstruct_array{cond_select},example_inds_cond_video(sorted_clust_ind(1:end),1),videodir,videotag)
animate_markers_aligned_fullmovie(mocapstruct_array{cond_select},example_inds_cond_video{sorted_clust_ind(130),1}(1:10:end));




figure(99999)
hist(number_cond(:,1)./sum(number_cond,2),100)

animate_markers_aligned_fullmovie(mocapstruct_caff,example_inds_cond_video{64,1}(1:10:end));
animate_markers_aligned_fullmovie(mocapstruct_prel,example_inds_cond_video{53,2});



% get_cand_frames()

%Need to refactor this
colors_plot = hsv(80);
colors_plot = colorcube(80)./2+lines(80)./2;
subset_of_points_to_plot2 = 1:100:100*6000;

legendnames = cell(1,0);

for mm =1:max(candidate_frames)
    [ framesubset_cand,framesubset_tsne,~] = intersect(subset_of_points_to_plot_embed,find(candidate_frames == mm));
    
    if numel(framesubset_cand)
        
        figure(572)
        plot(0, 0,'o','MarkerEdgeColor','none','MarkerSize',8,'MarkerFaceColor',colors_plot(mm,:))
        hold on
        legendnames{1,size(legendnames,2)+1} = fieldnames_beh{mm+1};
    end
end

for mm =1:max((candidate_frames))
    [ framesubset_cand,framesubset_tsne,~] = intersect(subset_of_points_to_plot_embed ,find(candidate_frames == mm));
    
    if numel(framesubset_cand)
        
        figure(572)
        plot( zValues(framesubset_tsne,1), zValues(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
        
        
    end
end

figure(559)
legend(legendnames)

figure(572)
legend(legendnames)


%
%             mappedX_eucold = tsne(candidate_features(subset_of_points_to_plot_tsne,:));
% mappedX_KL =  run_tSne(candidate_features(subset_of_points_to_plot_tsne,:),parameters,[]);
%



% sample points from the inverse prob. of the tsne maps

%