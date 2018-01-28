%% function to coembed two datasets
%s
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';

%% PARMAETERS
%% set parameters
%V8
ratname = 'JDM25';
conditionnumbers = [1,10];



ratname = 'Vicon8';
conditionname = 'Vicon8_caff';
conditionnumber = 6;

conditionname2 = 'Vicon8_prelesion';
conditionnumber2 = 1;

conditionname3 = 'Vicon8_amph';
conditionnumber3 = 7;
%% Load files

descriptor_struct = get_mocap_files_table(5,ratname);
 [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
 [mocapstruct_caff] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,1,[],mocapvideodirectory,1);

 descriptor_struct = get_mocap_files_table(conditionnumbers(2),ratname);
 [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
 [mocapstruct_prel] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,0,[],mocapvideodirectory,0);
 
  descriptor_struct = get_mocap_files_table(conditionnumber3,ratname);
 [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
 [mocapstruct_amph] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,0,[],mocapvideodirectory,0);
 
 
  ML_features_caff = get_supervised_features(mocapstruct_caff,mocapstruct_caff.modular_cluster_properties.clustering_inds_agg{2},2,ratname,'jdm25_caff',0);
ML_features_prel = get_supervised_features(mocapstruct_prel,mocapstruct_prel.modular_cluster_properties.clustering_inds_agg{2},2,ratname,'jdm25_prel',0);
ML_features_amph = get_supervised_features(mocapstruct_amph,mocapstruct_amph.modular_cluster_properties.clustering_inds_agg{2},2,ratname,'Vicon8_amph',0);

frames_with_goodtracking_caff = intersect(mocapstruct_caff.modular_cluster_properties.clustering_inds_agg{2},mocapstruct_caff.modular_cluster_properties.clipped_index{2});
frames_with_goodtracking_prel = intersect(mocapstruct_prel.modular_cluster_properties.clustering_inds_agg{2},mocapstruct_prel.modular_cluster_properties.clipped_index{2});
frames_with_goodtracking_amph = intersect(mocapstruct_amph.modular_cluster_properties.clustering_inds_agg{2},mocapstruct_amph.modular_cluster_properties.clipped_index{2});

[frames_with_goodtracking_moving_caff,ml_subset_moving_caff ] = intersect(mocapstruct_caff.move_frames,frames_with_goodtracking_caff);
[frames_with_goodtracking_moving_prel,ml_subset_moving_prel ] = intersect(mocapstruct_prel.move_frames,frames_with_goodtracking_prel);

%% Load annotations for each (Optional)

%%  run tsne on the individual files
  subset_of_points_to_plot_tsne = 1:100:100*20000;
    subset_of_points_to_plot_embed = 1:10:100*24000;

  
  tsne_features = {'spectrogram_pcs_wl_head_angle','spectrogram_pcs_wl_trunk_angle',...
      'ja_dyadic_spectrograms','appearance_features_agg_score_whitened'};
  num_feat = [10,10,10,6];
  
  candidate_features_caff_em = [];
    candidate_features_prel_em = [];
        candidate_features_amph_em = [];

    
 candidate_features_caff = [];
    candidate_features_prel = [];
        candidate_features_amph = [];

for mm = 1:numel(tsne_features)
      
      %      candidate_features_caff = cat(2,candidate_features_caff ,ML_features_caff.(tsne_features{mm})((subset_of_points_to_plot_tsne),1:num_feat(mm)));
   candidate_features_amph = cat(2,candidate_features_amph ,ML_features_amph.(tsne_features{mm})(subset_of_points_to_plot_tsne,1:num_feat(mm)));
           candidate_features_caff = cat(2,candidate_features_caff ,ML_features_caff.(tsne_features{mm})((subset_of_points_to_plot_tsne),1:num_feat(mm)));
     %candidate_features_prel = cat(2,candidate_features_prel ,ML_features_prel.(tsne_features{mm})(ml_subset_moving_prel(subset_of_points_to_plot_tsne),1:num_feat(mm)));
  
     candidate_features_caff_em = cat(2,candidate_features_caff_em ,ML_features_caff.(tsne_features{mm})((subset_of_points_to_plot_embed),1:num_feat(mm)));
    % candidate_features_prel_em = cat(2,candidate_features_prel_em ,ML_features_prel.(tsne_features{mm})(ml_subset_moving_prel(subset_of_points_to_plot_embed),1:num_feat(mm)));
      % candidate_features_caff_em = cat(2,candidate_features_caff_em ,ML_features_caff.(tsne_features{mm})((subset_of_points_to_plot_tsne),1:num_feat(mm)));
      candidate_features_amph_em = cat(2,candidate_features_amph_em ,ML_features_amph.(tsne_features{mm})(subset_of_points_to_plot_embed,1:num_feat(mm)));

end

parameters = [];
 %   mappedX_eucnew_caff =run_tSne(candidate_features_caff(subset_of_points_to_plot_tsne,:),parameters,'true');
  %  mappedX_eucnew_prel =run_tSne(candidate_features_caff(subset_of_points_to_plot_tsne,:),parameters,'true');

    %% smooth the tsne maps to sample from the watershed regions probabilistically
    numtosample = 10000;
%caff
% [xx,yy,density_caff] = findPointDensity(mappedX_eucnew_caff,2,[2001 2001],[-150 150]);
% reshaped_density_caff = reshape(density_caff,1,[]);
% sampled_caff = randsample(1:numel(reshaped_density_caff),numtosample,'true',reshaped_density_caff);
% 
% %prel
% [xx,yy,density_prel] = findPointDensity(mappedX_eucnew_prel,2,[2001 2001],[-150 150]);
% reshaped_density_prel = reshape(density_prel,1,[]);
% sampled_prel = randsample(1:numel(reshaped_density_prel),numtosample,'true',reshaped_density_prel);
inds_caff = randsample(1:numel(subset_of_points_to_plot_tsne),numtosample);
inds_prel = randsample(1:numel(subset_of_points_to_plot_tsne),numtosample);
inds_amph = randsample(1:numel(subset_of_points_to_plot_tsne),numtosample);

%%
jt_features = real(cat(1,candidate_features_caff(inds_caff,:),candidate_features_amph(inds_amph,:)));%candidate_features_amph(inds_amph,:)));%);%candidate_features_amph(inds_amph,:)); %,candidate_features_prel(inds_prel,:)
    mappedX_eucnew_joint =run_tSne(jt_features,parameters,'true');
  
%candidate_features_prel_em,

[zValues,outputStatistics] = ...
    findEmbeddings_precompfeatures(real(cat(1,candidate_features_caff_em,candidate_features_amph_em)),...
      jt_features,...
      mappedX_eucnew_joint,parameters,'true');
    
 
  
    
h=figure(605)
plot(mappedX_eucnew_joint(1:numel(inds_caff),1),mappedX_eucnew_joint(1:numel(inds_caff),2),'b+')
   hold on
   plot(mappedX_eucnew_joint(numel(inds_caff)+1:(numel(inds_caff)+numel(inds_amph)),1),...
       mappedX_eucnew_joint(numel(inds_caff)+1:(numel(inds_caff)+numel(inds_amph)),2),'r+')
   
   plot(mappedX_eucnew_joint((1+numel(inds_caff)+numel(inds_amph)):end,1),...
       mappedX_eucnew_joint(((numel(inds_caff)+numel(inds_amph)+1):end),2),'g+')
hold off

  
h= figure(605)
   hold on
plot(zValues(1:numel(subset_of_points_to_plot_embed),1),zValues(1:numel(subset_of_points_to_plot_embed),2),'bo','MarkerSize',1)
   plot(zValues(numel(subset_of_points_to_plot_embed)+1:2*numel(subset_of_points_to_plot_embed),1),...
       zValues(numel(subset_of_points_to_plot_embed)+1:2*numel(subset_of_points_to_plot_embed),2),'ro','MarkerSize',1)
     plot(zValues(2*numel(subset_of_points_to_plot_embed)+1:end,1),...
       zValues(2*numel(subset_of_points_to_plot_embed)+1:end,2),'go','MarkerSize',1)
hold off


%% track switching

[x, y] = getline(h);
  IN = inpolygon(zValues(:,1),zValues(:,2),...
        x ,y);
  %  good_track_switch = {frames_with_goodtracking_caff,frames_with_goodtracking_prel ,frames_with_goodtracking_amph};
   % cond_identifier = cat(1,ones(numel(subset_of_points_to_plot_embed),1),2*ones(numel(subset_of_points_to_plot_embed),1),...
   %     3*ones(numel(subset_of_points_to_plot_embed),1));
   %  cond_identifier = cat(1,ones(numel(subset_of_points_to_plot_embed),1),...
    %    3*ones(numel(subset_of_points_to_plot_embed),1));
            good_track_switch = {frames_with_goodtracking_caff,[],frames_with_goodtracking_amph};

      %  good_track_switch = {frames_with_goodtracking_moving_caff,frames_with_goodtracking_moving_prel };
 cond_identifier = cat(1,ones(numel(subset_of_points_to_plot_embed),1),...
        3*ones(numel(subset_of_points_to_plot_embed),1));
    
    
    beh_list = cell(1,3);
    for kk = [1 3]
    [~,ia] = intersect(find(cond_identifier==kk),find(IN));
    ia = subset_of_points_to_plot_embed(ia);
     beh_list{kk} = sort(reshape(unique(rectify_inds(bsxfun(@plus,good_track_switch{kk}(ia)',-10:10),max(good_track_switch{kk}))),1,[]),'ascend');    
    end
    close 370
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


    



%% look at embedded space velocity

figure(100)
plot(zValues)

zvel = sqrt(sum(abs(diff(zValues)).^2,2));
 zvel  = conv( zvel ,ones(1,5)./5,'same');
 
  figure(97)
    subplot(2,1,1)

  plot(zvel)
  subplot(2,1,2)
  
  hist(zvel,logspace(-4,4,1000))
  set(gca,'XScale','log')
  
  badvals = find(zvel>10);
    goodvals = find(zvel<10);

  figure(101)
  subplot(1,3,1)
  plot(zValues(badvals(1:end),1),zValues(badvals(1:end),2),'bo','MarkerSize',1)
  hold on
    plot(zValues(goodvals(1:end),1),zValues(goodvals(1:end),2),'ro','MarkerSize',1)
hold off
 subplot(1,3,2)
  plot(zValues(badvals(1:end),1),zValues(badvals(1:end),2),'bo','MarkerSize',1)
 subplot(1,3,3)
    plot(zValues(goodvals(1:end),1),zValues(goodvals(1:end),2),'ro','MarkerSize',1)

    

%% density, watershed, images
num_conditions = 3;
density_maps = cell(1,num_conditions);
[xx,yy,density_maps{1}] = findPointDensity(zValues(find(cond_identifier==1),:),2,[1001 1001],[-150 150]);
[xx,yy,density_maps{2}] = findPointDensity(zValues(intersect(goodvals,find(cond_identifier==2)),:),2,[1001 1001],[-150 150]);
[xx,yy,density_maps{3}] = findPointDensity(zValues(find(cond_identifier==3),:),2,[1001 1001],[-150 150]);

[xx,yy,density_jt] = findPointDensity(zValues(:,:),2,[1001 1001],[-150 150]);


  figure(480)
  subplot(1,4,1)
  imagesc(flipud(density_maps{1}))
    subplot(1,4,2)
  imagesc(flipud(density_maps{2}))
  subplot(1,4,3)
  imagesc(flipud(density_maps{3}))
   subplot(1,4,4)
  imagesc(flipud(density_jt))


  %% now threshold the density maps and apply the watershed
  figure(482)
 % [x,n] = hist(reshape(xcorr2(density_maps{3}),1,[]),1000);
  loglog(n,x)
  figure(582)
  [x,n] = hist(reshape((density_maps{3}),1,[]),logspace(-5,-3,200));
    [x2,n2] = hist(reshape((density_maps{1}),1,[]),logspace(-5,-3,200));

 loglog(n,x)
  hold on
  loglog(n2,x2,'r')
  hold off
  
  density_thresholds{1} = 3*10^(-5);
    density_thresholds{3} = 3*10^(-5);

  
  %figure(490)
  
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
imagesc(flipud(density_watersheds{3}))

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
 cond_select = 3;
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
   annotation_vec{ll} = zeros(1, numel(good_track_switch{ll}));
end

tsnehere = zValues;

subset_of_points_to_plot =1:size(tsnehere,1);
condition_labels = cat(1,ones(1,numel(subset_of_points_to_plot_embed)),3*ones(1,numel(subset_of_points_to_plot_embed)));

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

for ll = 1:3
[~,example_inds_cond{jj,ll}] = intersect(find(cond_identifier==ll),example_inds{jj});
annotation_vec{ll}(subset_of_points_to_plot_embed(example_inds_cond{jj,ll}))=jj;
end
number_cond(jj,:) = [numel(example_inds_cond{jj,1}) numel(example_inds_cond{jj,2}) numel(example_inds_cond{jj,3}) jj];




if numel(example_inds{jj})
temp = rectify_inds(unique(bsxfun(@plus,example_inds{jj},-10:10)),max(subset_of_points_to_plot_embed));
test = zeros(1,max(subset_of_points_to_plot_embed));
test(temp) = 1;
characteristics = bwconncomp(test);
cluster_num(jj) = characteristics.NumObjects;
cluster_mean_length(:,jj) = [characteristics.NumObjects mean(cellfun(@numel,characteristics.PixelIdxList)) ...
    median(cellfun(@numel,characteristics.PixelIdxList)) std(cellfun(@numel,characteristics.PixelIdxList))];

for kk =1:3
temp = (unique(bsxfun(@plus,good_track_switch{kk}(subset_of_points_to_plot_embed(example_inds_cond{jj,kk}))',-10:10)));
temp(temp<1) = 1;
        example_inds_cond_video{jj,kk} = temp;
end 
end
end
good_clusters = find(cellfun(@numel,example_inds_cond(:,3))>100)';
%sum()
make_dotplot(annotation_vec{3},annotation_labels,ML_features_amph,good_clusters)


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
      allframes =  cat(2,allframes,mocapstruct_amph.markers_aligned_preproc.(mocapstruct_amph.markernames{jj})(example_inds_cond_video{ll,3},:));
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
       dendrogram(Z,0)
       
figure(1111)
imagesc(pose_dist(sorted_clust_ind,sorted_clust_ind))
caxis([0 50])

L = density_watersheds{cond_select};
s = regionprops(L, 'Centroid');
hold on
ind = 0;
for k = sorted_clust_ind'
    ind = ind+1;
          L(density_cc{cond_select}.PixelIdxList{k}) = ind ;
end

      figure(486)
   imagesc(L)
   ind = 0;

for k = sorted_clust_ind'
    ind = ind+1;
    c = s(k).Centroid;
    text(c(1), c(2), num2str(ind), ... %sprintf('%d', integer(ind)),
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','Color','white');
end
hold off
   textColor = 'white';
   


       
plot_timelapse_mosaic(mocapstruct_caff,example_inds_cond_video(sorted_clust_ind(1:end),1))
  animate_markers_aligned_fullmovie(mocapstruct_amph,example_inds_cond_video{sorted_clust_ind(111),3}(1:10:end));




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