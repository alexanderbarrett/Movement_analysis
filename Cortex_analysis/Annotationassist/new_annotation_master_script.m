%% set this to the server address
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';

%% PARMAETERS
%for the ML features
ratname = 'Vicon8';
conditionnumbers = 6;
annotation_number = 2;

%JDM25
% ratname = 'JDM25';
%      conditionnumbers = 5; %1 is pre lesion, 10 is post lesion
%      
    
     %JDM33
     %ratname = 'JDM33';
     %conditionnumber = [1,5,10]; %1 is pre lesion, 10 is post lesion
     

% annotation_number = 4; % look in the file below to check which file to load in

[fullposture_annotation_struct,fulloutput_annotation_struct,subset_of_frames_annotated] = load_handannotation_files(annotation_number,mocapmasterdirectory);

%% get files and the machine learning features
descriptor_struct = get_mocap_files_table(conditionnumbers(1),ratname,mocapmasterdirectory);
 [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
 [mocapstruct_all] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,0,[],mocapvideodirectory,1);
ML_features = get_supervised_features(mocapstruct_all,mocapstruct_all.modular_cluster_properties.clustering_inds_agg{2},2,ratname,descriptor_struct.Nametag,0,0);


%% visualize a chosen behavior (optional)
%animate_markers_aligned_fullmovie(mocapstruct_all,fulloutput_annotation_struct.WetDogShake(1:10:end)')

%% Convert the struct to a categorical output for the random forest, ie each
% behavior is a number
subset_of_frames_annotated_full = 1:max(subset_of_frames_annotated);
[outputvector,fieldnames_beh,observed_behaviors] = generate_categorical_output(fulloutput_annotation_struct,max(subset_of_frames_annotated));
fieldnames_beh = cat(1,'null',fieldnames_beh); %there is a 0 for unannotated in the categorical vector

%% remove unannotated frames
[~,ind_unannotated] = intersect(subset_of_frames_annotated_full,find(outputvector == 0));
subset_of_frames_annotated_full(ind_unannotated)= [];

% get the outputvector- intersect the frames that have good tracking (
[ind_for_outputvec,annotated_frames_goodtracking,~] = intersect(ML_features.framelist_true,subset_of_frames_annotated_full);
annotation_output_goodtracking = outputvector(ind_for_outputvec);
% get the subset of overall frames with good tracking
%[~,annotated_frames_goodtracking] = intersect(frames_with_goodtracking,subset_of_frames_annotated);


%% RUN THE ML ALGORITHM TO FIND SIMILAR FRAMES
fprintf('searching for similar frames \n')
frames_to_annotate = 1:max(subset_of_frames_annotated);
clear candidate_frames
clear score
clear Mdl_test
[candidate_frames,score,Mdl_test,annotated_instances] = findsimilarframes_mlfeatures_multifeatures(ML_features,annotation_output_goodtracking,annotated_frames_goodtracking,frames_to_annotate);
% [candidate_frames,score,Mdl_test] = findsimilarframes_mlfeatures_multifeatures(ML_features,labels,frames,frames_to_annotate);

oob_err = oobError(Mdl_test,'mode','ensemble') ;

[all_behaviors,cfmat,tp_rate] = beh_analysis(annotated_instances,candidate_frames,fieldnames_beh,observed_behaviors);

% To display all of the behaviors in command window at once
behs = fieldnames(all_behaviors);
for i = 1:numel(behs)  
  if all_behaviors.(behs{i}).num_annotated ~= 0   
      fprintf(behs{i})
      fprintf('\n')
      display(all_behaviors.(behs{i})) 
  end
end

behs = fieldnames(all_behaviors);
ann = 0;
pred = 0; 
corr = 0; 
for i = 1:numel(behs)  
  ann = ann + all_behaviors.(behs{i}).num_annotated ; 
  pred = pred + all_behaviors.(behs{i}).num_predicted ; 
  corr = corr + all_behaviors.(behs{i}).num_correct_pred ;     
end


%% VISUALIZE AN EXAMPLE BEHAVIOR
ML_framelist = ML_features.framelist_true;

chosenbehavior = 'FaceWipes';
animate_markers_aligned_fullmovie(mocapstruct_all,fulloutput_annotation_struct.(chosenbehavior)(1:10:end)')

candidate_frames_vec = (find(candidate_frames==find(strcmp(fieldnames_beh ,chosenbehavior)==1)-1));
candidate_frame_predicted = setxor(ML_framelist(candidate_frames_vec),...
    intersect(ML_framelist,fulloutput_annotation_struct.(chosenbehavior)));
%predicted_annotatedframes= candidate_frames_vec(candidate_frames_vec<80000);
%predicted_unannotatedframes= candidate_frames_vec(candidate_frames_vec>80000);

%visualize the ground truth and predicted data
animate_markers_aligned_fullmovie(mocapstruct_all,(candidate_frame_predicted(1:10:end))')
%animate_markers_aligned_fullmovie(mocapstruct_all,frames_with_goodtracking(predicted_annotatedframes(1:10:end))')



%% look at behavior characteristics

%get only the types of behavior that are found
unique_predictions = unique(candidate_frames(~isnan(candidate_frames)));
%these are the fieldnames of these behaviors
fieldnames_beh(unique_predictions+1)



        make_dotplot(    candidate_frames,fieldnames_beh,ML_features,unique_predictions)
                make_dotplot(    outputvector,fieldnames_beh,ML_features,unique_predictions)
[agg_dprimes,agg_names] =  get_feature_dprimes(ML_features,outputvector,unique_predictions,fieldnames_beh);
[transout,statefreq] = compute_transition_matrix(candidate_frames,fieldnames_beh);


%% add other functions to compare real and found

%visualize the behavioral segregation
%% plot on the tsne
%Need to refactor this

subset_of_points_to_plot2 = 1:100:100*6000;
subset_of_points_to_plot = 1:200:200*6000;
subset_of_points_to_plot= 1:100:100*15000;






%tsnehere =  mapped_jt_alldynamics_angles;
tsnehere = map;
[xx,yy,density] = findPointDensity( tsnehere,1,[2001 2001],[-50 50]);
figure(573)
subplot(1,2,1)
imagesc(flipud(density))
colorbar
subplot(1,2,2)
inv_density = 1./density;
inv_density(inv_density>10^6) = 10^10;
L=watershed(inv_density);
imagesc(flipud(L))
conn_watershed = bwconncomp(L);
   stats = regionprops(conn_watershed,'convexhull');

% % 
% thresh = 1*10^(-4);
%  density(density<thresh) = 0;
%  imagesc(density)
% [conn_watershed] = bwconncomp(density);
%    stats = regionprops(conn_watershed,'convexhull');


example_inds = cell(1,conn_watershed.NumObjects);
example_inds_exp = cell(1,conn_watershed.NumObjects);
cluster_num = zeros(1,conn_watershed.NumObjects);
cluster_mean_length = zeros(4,conn_watershed.NumObjects);


for jj = 2:conn_watershed.NumObjects
   [example_inds_x,example_inds_y]  = ind2sub(conn_watershed.ImageSize,conn_watershed.PixelIdxList{jj} );
   verts_x = (round(stats(jj).ConvexHull(:,1)));
verts_x(verts_x<1) = 1;
   verts_x(verts_x>numel(xx)) =numel(xx);
   
    verts_y = (round(stats(jj).ConvexHull(:,2)));
   verts_y(verts_y<1) = 1;
   verts_y(verts_y>numel(yy)) =numel(yy);
   %% get the points in the watershed polygon
   IN = inpolygon(tsnehere(:,1),tsnehere(:,2),...
        xx(verts_x) ,yy(verts_y));
    
    %% get the bout characteristics
example_inds{jj} = subset_of_points_to_plot(find(IN));
temp = (unique(bsxfun(@plus,example_inds{jj}',-10:10)));
temp(temp<1) = 1;
temp(temp>max(subset_of_points_to_plot)) = max(subset_of_points_to_plot);
test = zeros(1,max(subset_of_points_to_plot));
test(temp) = 1;
characteristics = bwconncomp(test);
cluster_num(jj) = characteristics.NumObjects;
cluster_mean_length(:,jj) = [characteristics.NumObjects mean(cellfun(@numel,characteristics.PixelIdxList)) ...
    median(cellfun(@numel,characteristics.PixelIdxList)) std(cellfun(@numel,characteristics.PixelIdxList))];


temp = (unique(bsxfun(@plus,example_inds{jj}',-20:20)));
temp(temp<1) = 1;
temp(temp>numel(frames_with_goodtracking)) = numel(frames_with_goodtracking);
temp = frames_with_goodtracking(temp);

example_inds_exp{jj} = temp(1:10:end);
end

plot_pose_mosaic(mocapstruct_all,example_inds_exp)
M = plot_multi_clusters(mocapstruct_all,example_inds_exp,50);

indshere = frames_with_goodtracking(unique(bsxfun(@plus,example_inds{20}',-50:50)));
animate_markers_aligned_fullmovie(mocapstruct_all,example_inds_exp{180} );


  v = VideoWriter(strcat('aggregate_movie_4'),'MPEG-4');
    open(v)
    writeVideo(v, M)
    close(v)





mappedX_joint =  tsne(cat(2,ML_features.pose_score(subset_of_points_to_plot,1:10),ML_features.appearance_features_agg_score_whitened(subset_of_points_to_plot,1:6),...
    ML_features.spectrogram_pcs_trunk_angle(subset_of_points_to_plot,1:15),ML_features.spectrogram_pcs_head_angle(subset_of_points_to_plot,1:15)));

mappedX = tsne(cat(2,ML_features.pose_score(subset_of_points_to_plot,1:10),ML_features.appearance_features_agg_score_whitened(subset_of_points_to_plot,1:6)));

mappedX_japose = tsne(ML_features.ja_dyadic_spectrograms(subset_of_points_to_plot,1:10));

 mapped_trunk = tsne(cat(2,ML_features.spectrogram_pcs_wl_trunk_angle(subset_of_points_to_plot2,1:10)));
   
   mapped_jt_dynamics = tsne(cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot,1:15),ML_features.spectrogram_pcs_wl_trunk_angle(subset_of_points_to_plot,1:15)));
   
   % mapped_all_dynamics = tsne(cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot,1:15),...
   %     ML_features.spectrogram_pcs_wl_trunk_angle(subset_of_points_to_plot,1:15)));
  
 mapped_head = tsne(cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot2,1:15)));
 
      mappedX_dyn_angle_t =  tsne(cat(2,ML_features.spectrogram_pcs_trunk_angle(subset_of_points_to_plot2,1:15)));
    
        mapped_jt_alldynamics = tsne(cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot2,1:15),...
          ML_features.spectrogram_pcs_wl_trunk_angle(subset_of_points_to_plot2,1:15),...
          ML_features.spectrogram_pcs_head_angle(subset_of_points_to_plot2,1:15),...
          ML_features.spectrogram_pcs_trunk_angle(subset_of_points_to_plot2,1:15)));
      
        % not great
      mapped_jt_alldynamics_angles_pose = tsne(cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot,1:15),...
          ML_features.spectrogram_pcs_wl_trunk_angle(subset_of_points_to_plot,1:15),...
          ML_features.ja_dyadic_spectrograms(subset_of_points_to_plot,1:10),...
      ML_features.pose_score(subset_of_points_to_plot,1:10),...
      ML_features.appearance_features_agg_score_whitened(subset_of_points_to_plot,1:6),...
          ML_features.spectrogram_pcs_head_angle(subset_of_points_to_plot,1:15),...
          ML_features.spectrogram_pcs_trunk_angle(subset_of_points_to_plot,1:15)));
      



  % okay
   mapped_jt_alldynamics_angles_only = tsne(cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot,1:15),...
          ML_features.spectrogram_pcs_wl_trunk_angle(subset_of_points_to_plot,1:15),...
          ML_features.ja_dyadic_spectrograms(subset_of_points_to_plot,1:10),...
          ML_features.appearance_features_agg_score_whitened(subset_of_points_to_plot,1:6) )...      
          );
   
mappedX_dyn_head=  tsne(cat(2,ML_features.spectrogram_pcs_head_angle(subset_of_points_to_plot,1:15),ML_features.spectrogram_pcs_trunk_angle(subset_of_points_to_plot,1:15)));

    
    
      %% try the embedding
      subset_of_points_to_plot_tsne = 1:100:100*10000;
subset_of_points_to_plot = 1:30:300*7000;
      subset_of_points_to_plot_emreal = 1:3:300*7000;

      % tsne on a subset
      chosenbehavior = 'Walk';
candidate_frames_vec = (find(candidate_frames==find(strcmp(fieldnames_beh ,chosenbehavior)==1)-1));
   subset_of_points_to_plot =candidate_frames_vec(1:10:end);

   
           % will try - worked well before
           feat_vec_pose = cat(2,...
          ML_features.ja_dyadic_spectrograms(:,1:10),...
      ML_features.pose_score(:,1:10),...
      ML_features.appearance_features_agg_score_whitened(:,1:6));
           
     feat_vec_dyn = cat(2,cat(1,ML_features.spectrogram_pcs_wl_head_angle(:,1:15),zeros(1,15)),...
          cat(1,ML_features.spectrogram_pcs_wl_trunk_angle(:,1:15),zeros(1,15)),...
                    cat(1,ML_features.ja_eig_wlpcs(:,1:15),zeros(1,15)),...
      ML_features.absolute_velocity_trunk_abs_100(:),ML_features.absolute_std_velocity_trunk_abs_100(:),...
        ML_features.rel_velocity_hipR_abs_100(:),ML_features.rel_velocity_hipL_abs_100(:),...
        ML_features.rel_std_velocity_hipR_abs_100(:),ML_features.rel_std_velocity_hipL_abs_100(:),...
   ML_features.rel_velocity_head_abs_100(:),ML_features.rel_std_velocity_head_abs_100(:),...
   ML_features.rel_velocity_trunk_abs_100(:),ML_features.rel_std_velocity_trunk_abs_100(:),...
    ML_features.absolute_velocity_trunk_abs_300(:),ML_features.absolute_std_velocity_trunk_abs_300(:),...
        ML_features.rel_velocity_hipR_abs_300(:),ML_features.rel_velocity_hipL_abs_300(:),...
        ML_features.rel_std_velocity_hipR_abs_300(:),ML_features.rel_std_velocity_hipL_abs_300(:),...
   ML_features.rel_velocity_head_abs_300(:),ML_features.rel_std_velocity_head_abs_300(:),...
   ML_features.rel_velocity_trunk_abs_300(:),ML_features.rel_std_velocity_trunk_abs_300(:),...
      ML_features.absolute_velocity_trunk_abs_33(:),ML_features.absolute_std_velocity_trunk_abs_33(:),...
        ML_features.rel_velocity_hipR_abs_33(:),ML_features.rel_velocity_hipL_abs_33(:),...
        ML_features.rel_std_velocity_hipR_abs_33(:),ML_features.rel_std_velocity_hipL_abs_33(:),...
   ML_features.rel_velocity_head_abs_33(:),ML_features.rel_std_velocity_head_abs_33(:),...
   ML_features.rel_velocity_trunk_abs_33(:),ML_features.rel_std_velocity_trunk_abs_33(:),...
   ML_features.ja_velocitytrunk_angle_difforder_100, ML_features.ja_velocityhead_angle_difforder_100,...
     ML_features.ja_stdtrunk_angle_difforder_100, ML_features.stdhead_angle_difforder_100 );
  
  
           feat_vec = cat(2,cat(1,ML_features.spectrogram_pcs_wl_head_angle(:,1:15),zeros(1,15)),...
          cat(1,ML_features.spectrogram_pcs_wl_trunk_angle(:,1:15),zeros(1,15)),...
          ML_features.ja_dyadic_spectrograms(:,1:10),...
      ML_features.pose_score(:,1:10),...
      ML_features.appearance_features_agg_score_whitened(:,1:6),...
      ML_features.absolute_velocity_trunk_abs_100(:),ML_features.absolute_std_velocity_trunk_abs_100(:),...
        ML_features.rel_velocity_hipR_abs_100(:),ML_features.rel_velocity_hipL_abs_100(:),...
        ML_features.rel_std_velocity_hipR_abs_100(:),ML_features.rel_std_velocity_hipL_abs_100(:),...
   ML_features.rel_velocity_head_abs_100(:),ML_features.rel_std_velocity_head_abs_100(:),...
   ML_features.rel_velocity_trunk_abs_100(:),ML_features.rel_std_velocity_trunk_abs_100(:),...
    ML_features.absolute_velocity_trunk_abs_300(:),ML_features.absolute_std_velocity_trunk_abs_300(:),...
        ML_features.rel_velocity_hipR_abs_300(:),ML_features.rel_velocity_hipL_abs_300(:),...
        ML_features.rel_std_velocity_hipR_abs_300(:),ML_features.rel_std_velocity_hipL_abs_300(:),...
   ML_features.rel_velocity_head_abs_300(:),ML_features.rel_std_velocity_head_abs_300(:),...
   ML_features.rel_velocity_trunk_abs_300(:),ML_features.rel_std_velocity_trunk_abs_300(:),...
      ML_features.absolute_velocity_trunk_abs_33(:),ML_features.absolute_std_velocity_trunk_abs_33(:),...
        ML_features.rel_velocity_hipR_abs_33(:),ML_features.rel_velocity_hipL_abs_33(:),...
        ML_features.rel_std_velocity_hipR_abs_33(:),ML_features.rel_std_velocity_hipL_abs_33(:),...
   ML_features.rel_velocity_head_abs_33(:),ML_features.rel_std_velocity_head_abs_33(:),...
   ML_features.rel_velocity_trunk_abs_33(:),ML_features.rel_std_velocity_trunk_abs_33(:));%,...
 %   ML_features.ja_velocitytrunk_angle_difforder_100, ML_features.ja_velocityhead_angle_difforder_100,...
  %   ML_features.ja_stdtrunk_angle_difforder_100, ML_features.ja_stdhead_angle_difforder_100 ,...
     %          cat(1,ML_features.ja_eig_wlpcs(:,1:15),zeros(1,15)));

%feat_vec = feat_vec_dyn;
feat_vec = bsxfun(@rdivide,bsxfun(@minus,feat_vec,mean(feat_vec,1)),nanstd(feat_vec,[],1));

numDims = 2; pcaDims = 30; perplexity = 30; theta = .5; alg = 'svd';
  tic
       map =  fast_tsne(feat_vec(subset_of_points_to_plot,:), numDims, pcaDims, perplexity, theta, alg);
        toc



        mapped_jt_alldynamics_angles_pose_noxyz = run_tSne(feat_vec(subset_of_points_to_plot,:),[],'true');

      
        
        figure(909)
        hold on
plot(map(:,1),map(:,2),'+b')
    
        
            mappedX_eucold = tsne(candidate_features(subset_of_points_to_plot_tsne,:));
            
mappedX_KL =  run_tSne(candidate_features(subset_of_points_to_plot_tsne,:),parameters,[]);
    mappedX_eucnew =run_tSne(candidate_features(subset_of_points_to_plot_tsne,:),parameters,'true');
    

[zValues,outputStatistics] = ...
    findEmbeddings_precompfeatures(feat_vec(subset_of_points_to_plot_emreal,:),...
      feat_vec(subset_of_points_to_plot,:),...
      map,[],'true');

  subset_of_points_to_plot_emreal = subset_of_points_to_plot_emreal(  find(outputStatistics.inConvHull==1));
zValuesin =zValues(  find(outputStatistics.inConvHull==1),:);
  
figure(999)
plot(zValues)
hold on
plot(outputStatistics.zGuesses,'r')
hold off
  
  zvel = sqrt(sum(abs(diff(zValues)).^2,2));
 zvel  = conv( zvel ,ones(1,3)./3,'same');
 
  figure(97)
    subplot(2,1,1)

  plot(zvel)
  subplot(2,1,2)
  
  hist(zvel,logspace(-4,4,1000))
  set(gca,'XScale','log')
  
  badvals = find(zvel>10);
    goodvals = find(zvel<1);

  figure(101)
  plot(zValues(badvals,1),zValues(badvals,2),'b+')
  hold on
    plot(zValues(goodvals,1),zValues(goodvals,2),'r+')
hold off
  
  
  framesanim = frames_with_goodtracking(subset_of_points_to_plot_emreal(badvals));
framesanim = unique(bsxfun(@plus,framesanim',-20:20));
  animate_markers_aligned_fullmovie(mocapstruct_all,framesanim(1:10:end) );
  
figure(99)
plot3(ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot,1),...
    ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot,2),...
    ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot,3),'+')

%get a tsne if desired
%mappedX=  tsne(cat(2,ML_features.spectrogram_pcs_head_angle(subset_of_points_to_plot,1:10),...
 %   ML_features.spectrogram_pcs_trunk_angle(subset_of_points_to_plot,1:10)));

 colors_plot = parula(80);
%colors_plot = colorcube(80)./2+lines(80)./2;
 
colors_plot =cbrewer('qual','Paired', numel(unique(candidate_frames)));

legendnames = cell(1,0);
beh_search =unique(candidate_frames);
for mm =beh_search'
       [ framesubset_cand,framesubset_tsne,~] = intersect(subset_of_points_to_plot,find(candidate_frames == mm));
       
   if numel(framesubset_cand)
figure(559)
plot3(0,0, 0,'o','MarkerEdgeColor','none','MarkerSize',8,'MarkerFaceColor',colors_plot(find(beh_search==mm),:))
hold on

figure(572)
plot(0, 0,'o','MarkerEdgeColor','none','MarkerSize',8,'MarkerFaceColor',colors_plot(find(beh_search==mm),:))
hold on
        legendnames{1,size(legendnames,2)+1} = fieldnames_beh{mm+1};
   end
end

for mm =beh_search'
   [ framesubset_cand,framesubset_tsne,~] = intersect(subset_of_points_to_plot,find(candidate_frames == mm));

   if numel(framesubset_cand)
figure(559)
plot3(ML_features.pose_score(framesubset_cand,1),ML_features.pose_score(framesubset_cand,2),...
   ML_features.pose_score(framesubset_cand,3),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(find(beh_search==mm),:))

% plot3(ML_features.ja_dyadic_spectrograms(framesubset_cand,1),ML_features.ja_dyadic_spectrograms(framesubset_cand,2),...
%     ML_features.ja_dyadic_spectrograms(framesubset_cand,3),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))

%plot3(ML_features.spectrogram_pcs_trunk_angle(framesubset_cand,1),ML_features.spectrogram_pcs_trunk_angle(framesubset_cand,2),...
 %   ML_features.spectrogram_pcs_trunk_angle(framesubset_cand,3),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))


%ML_features.ja_dyadic_spectrograms;

figure(572)
%plot(mappedX_dyn_angle(framesubset_tsne,1),mappedX_dyn_angle(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
%plot(mappedX_dyn_head(framesubset_tsne,1),mappedX_dyn_head(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
plot( map(framesubset_tsne,1), map(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(find(beh_search==mm),:))


   end
end

 figure(559)
legend(legendnames)

 h=figure(909)
legend(legendnames)


zValues =map;

[x, y] = getline(h);
  IN = inpolygon(zValues(:,1),zValues(:,2),x ,y);

animate_markers_aligned_fullmovie(mocapstruct_caff,unique(bsxfun(@plus,frames_with_goodtracking(subset_of_points_to_plot(find(IN)))',0)))


make_ethogram(outputvector,fieldnames_beh)
plot_behavior_examples(mocapstruct_all,fulloutput_annotation_struct.GeneralGroom,indivbouts_annotation_struct.GeneralGroom)

% make example PSDs etc. 
legend_tags = {'caff','bidls','bimc'};
cellstruct = cell(1,3);
cellstruct{1} = mocapstruct_caff;
cellstruct{2} = mocapstruct_concatenated;
cellstruct{3} = mocapstruct_post_bi;
timestruct{1} = fulloutput_caff;
timestruct{2} = fulloutput_DLS2;
timestruct{3} = fulloutput_mc;


for kk =1:3
multi_plot_marker_characteristics_timerange(cellstruct{kk},timestruct{kk},colorarray{kk},kk)
end
