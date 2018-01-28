
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';

%% PARMAETERS
%for the ML features
ratname = 'Vicon8';
conditionname = 'Vicon8_caff';
conditionnumber = 6;

conditionname = 'Vicon8_prelesion';
conditionnumber = 1;

% for the specific file
subset_of_frames_to_annotate = 1:2.5*10^6; %subset of well tracked frames to annotate
gapfill_number = 20;

%% can load in multiple filenames
%annotation_mocapname =
%'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Preprocessed\nolj_Recording_day7_caffeine1_nolj.mat';.
%% FIRST FILE
annotation_folder = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Preprocessed\';

annotation_filenames = {'nolj_Recording_day7_caffeine1_nolj_handannotation_JDM.mat'};
annotation_filenumber = {1};

%% SECOND GROUP OF FILES
annotation_folder = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Preprocessed\';
annotation_filenames = {'nolj_Recording_day7_caffeine1_nolj_handannotation_JDM.mat',...
    'nolj_Recording_day7_caffeine1_nolj_handannotation_AB_1to48000.mat',...
    'nolj_Recording_day7_caffeine2_nolj_handannotation_AB_1to23800.mat',...
    'nolj_Recording_day7_caffeine2_nolj_handannotation_AB_23800to53323.mat',...
    'nolj_Recording_day7_caffeine2_nolj_handannotation_AB_53323to64603.mat'};

annotation_filenumber = {1,1,2,2,2};
%% THIRD GROUP 
%annotation_folder = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170820\Preprocessed\';
%annotation_filenames = {'Recording_day5_overnight27_nolj_handannotation2_AB.mat'
annotation_filenumber = {23};

annotation_folder = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170817\Preprocessed\';

annotation_filenames = {'Recording_overnight_day28_nolj_handannotation_AB',...
    'Recording_overnight_day29_nolj_handannotation_AB',...
    'Recording_overnight_day211_nolj_handannotation_AB',...
    'Recording_overnight_day215_nolj_handannotation_AB',...
        'Recording_overnight_day216_nolj_handannotation_AB',...
        'Recording_overnight_day218_nolj_handannotation_AB',...
        'Recording_overnight_day220_nolj_handannotation_AB'...
    };
annotation_filenumber = {1,2,4,8,9,11,13};

annot_cell = cell(1,numel(annotation_filenames));
subset_of_frames_annotated = [];
agg_struct = struct([]);
agg_posture = struct([]);

for jj = 1:numel(annotation_filenames)%2:5
output1 = load(strcat(annotation_folder,annotation_filenames{jj}));
  C = struct2cell(output1.output.GlobalBehavior);
minval = min([C{:}]);
maxval = max([C{:}]);

for kk = fieldnames(output1.output.GlobalBehavior)'
    output1.output.GlobalBehavior.(kk{1}) = output1.output.GlobalBehavior.(kk{1})+(annotation_filenumber{jj}-1)*540000;
end
for kk = fieldnames(output1.output.Posture)'
        output1.output.Posture.(kk{1}) = output1.output.Posture.(kk{1})+(annotation_filenumber{jj}-1)*540000;
end
annot_cell{jj} = output1;


subset_of_frames_annotated = cat(2,subset_of_frames_annotated ,(annotation_filenumber{jj}-1)*540000+(minval:maxval)); %depends on the file(s) loaded
agg_struct = mergeStructs_JDM(agg_struct,output1.output.GlobalBehavior);
agg_posture = mergeStructs_JDM(agg_posture,output1.output.Posture);

end
subset_of_frames_annotated = unique(subset_of_frames_annotated);

[pose_struct_out,globalbehavior_struct_out] = posture_behavior_gapfill(agg_posture,agg_struct);
%outputstruct_annotation = load(annotation_filename);


%% get files and the machine learning features
descriptor_struct = get_mocap_files_table(conditionnumber,ratname);
 [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
 [mocapstruct_all] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,0,[],mocapvideodirectory,0);
ML_features = get_supervised_features(mocapstruct_all,mocapstruct_all.modular_cluster_properties.clustering_inds_agg{2},2,ratname,'Vicon8_caff',0);
%the frames are the intersection of those that are returned from clipping
%and in the cluster of markers
frames_with_goodtracking = intersect(mocapstruct_all.modular_cluster_properties.clustering_inds_agg{2},mocapstruct_all.modular_cluster_properties.clipped_index{2});

%% LOAD IN THE ANNOTATION
% fill the gaps in the structure
% in the future can write code to accomodate multiple annotation_files 
[fulloutput_annotation_struct,indivbouts_annotation_struct] = fillannotationgaps_struct(agg_struct ,gapfill_number);
[fullposture_annotation_struct,indivbouts_posture_struct] = fillannotationgaps_struct(agg_posture ,gapfill_number);

%visualize a chosen behavior (optional)
%animate_markers_aligned_fullmovie(mocapstruct_caff,fulloutput_caff_struct.(chosenbehavior)(1:10:end)')

% Convert the struct to a categorical output for the random forest
[outputvector,fieldnames_beh,observed_behaviors] = generate_categorical_output(fulloutput_annotation_struct,numel(subset_of_frames_annotated));
fieldnames_beh = cat(1,'null',fieldnames_beh); %there is a 0 for unannotated in the categorical vector

%% remove unannotated frames
[~,ind_unannotated] = intersect(subset_of_frames_annotated,find(outputvector == 0));
subset_of_frames_annotated(ind_unannotated)= [];

% get the outputvector
[ind_for_outputvec,annotated_frames_goodtracking,~] = intersect(frames_with_goodtracking,subset_of_frames_annotated);
annotation_output_goodtracking = outputvector(ind_for_outputvec);
% get the subset of overall frames with good tracking
%[~,annotated_frames_goodtracking] = intersect(frames_with_goodtracking,subset_of_frames_annotated);

%% RUN THE ML ALGORITHM TO FIND SIMILAR FRAMES
fprintf('searching for similar frames \n')
[candidate_frames,score,Mdl_test] = findsimilarframes_mlfeatures_multifeatures(ML_features,annotation_output_goodtracking,annotated_frames_goodtracking,subset_of_frames_to_annotate);
%get only the types of behavior that are found
unique_predictions = unique(candidate_frames(~isnan(candidate_frames)));
%these are the fieldnames of these behaviors
fieldnames_beh(unique_predictions+1)

%% VISUALIZE AN EXAMPLE BEHAVIOR
animate_markers_aligned_fullmovie(mocapstruct_all,fulloutput_annotation_struct.(chosenbehavior)(1:10:end)')
chosenbehavior = 'FaceWipes';

candidate_frames_vec = (find(candidate_frames==find(strcmp(fieldnames_beh ,chosenbehavior)==1)-1));
candidate_frame_predicted = setxor(frames_with_goodtracking(candidate_frames_vec),...
    intersect(frames_with_goodtracking,fulloutput_annotation_struct.(chosenbehavior)));
%predicted_annotatedframes= candidate_frames_vec(candidate_frames_vec<80000);
%predicted_unannotatedframes= candidate_frames_vec(candidate_frames_vec>80000);

%visualize the ground truth and predicted data
animate_markers_aligned_fullmovie(mocapstruct_all,(candidate_frame_predicted(1:10:end))')
%animate_markers_aligned_fullmovie(mocapstruct_all,frames_with_goodtracking(predicted_annotatedframes(1:10:end))')

%% add other functions to compare real and found

%visualize the behavioral segregation
%% plot on the tsne
%Need to refactor this
colors_plot = hsv(80);
colors_plot = colorcube(80)./2+lines(80)./2;
subset_of_points_to_plot2 = 1:100:100*6000;
subset_of_points_to_plot = 1:100:100*20000;
subset_of_points_to_plot= 1:100:100*15000;

%tsnehere =  mapped_jt_alldynamics_angles;
tsnehere = zValues;
[xx,yy,density] = findPointDensity( tsnehere,2,[2001 2001],[-150 150]);
figure(573)
subplot(1,2,1)
imagesc(density)
colorbar
subplot(1,2,2)
inv_density = 1./density;
inv_density(inv_density>10^6) = 10^10;
L=watershed(inv_density);
imagesc(L)
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
      
           % will try - worked well before
        mapped_jt_alldynamics_angles_pose_noxyz = tsne(cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot,1:15),...
          ML_features.spectrogram_pcs_wl_trunk_angle(subset_of_points_to_plot,1:15),...
          ML_features.ja_dyadic_spectrograms(subset_of_points_to_plot,1:10),...
      ML_features.pose_score(subset_of_points_to_plot,1:10),...
      ML_features.appearance_features_agg_score_whitened(subset_of_points_to_plot,1:6)...
        ));

    
  % okay
   mapped_jt_alldynamics_angles_only = tsne(cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot,1:15),...
          ML_features.spectrogram_pcs_wl_trunk_angle(subset_of_points_to_plot,1:15),...
          ML_features.ja_dyadic_spectrograms(subset_of_points_to_plot,1:10),...
          ML_features.appearance_features_agg_score_whitened(subset_of_points_to_plot,1:6) )...      
          );
   
mappedX_dyn_head=  tsne(cat(2,ML_features.spectrogram_pcs_head_angle(subset_of_points_to_plot,1:15),ML_features.spectrogram_pcs_trunk_angle(subset_of_points_to_plot,1:15)));

    
    
      %% try the embedding
      subset_of_points_to_plot_tsne = 1:100:100*10000;
subset_of_points_to_plot_em = 1:1:100*10000;
      subset_of_points_to_plot_emreal = 1:10:100*10000;

      candidate_features = cat(2,ML_features.spectrogram_pcs_wl_head_angle(subset_of_points_to_plot_em,1:15),...
          ML_features.spectrogram_pcs_wl_trunk_angle(subset_of_points_to_plot_em,1:15),...
          ML_features.ja_dyadic_spectrograms(subset_of_points_to_plot_em,1:10),...
      ML_features.pose_score(subset_of_points_to_plot_em,1:10),...
      ML_features.appearance_features_agg_score_whitened(subset_of_points_to_plot_em,1:6) );
      parameters = [];

            mappedX_eucold = tsne(candidate_features(subset_of_points_to_plot_tsne,:));
mappedX_KL =  run_tSne(candidate_features(subset_of_points_to_plot_tsne,:),parameters,[]);
    mappedX_eucnew =run_tSne(candidate_features(subset_of_points_to_plot_tsne,:),parameters,'true');
    

[zValues,outputStatistics] = ...
    findEmbeddings_precompfeatures(candidate_features(subset_of_points_to_plot_emreal,:),...
      candidate_features(subset_of_points_to_plot_tsne,:),...
      mappedX_eucnew,parameters,'true');

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

legendnames = cell(1,0);

for mm =1:max(candidate_frames)
       [ framesubset_cand,framesubset_tsne,~] = intersect(subset_of_points_to_plot_emreal,find(candidate_frames == mm));
       
   if numel(framesubset_cand)
figure(559)
plot3(0,0, 0,'o','MarkerEdgeColor','none','MarkerSize',8,'MarkerFaceColor',colors_plot(mm,:))
hold on

figure(572)
plot(0, 0,'o','MarkerEdgeColor','none','MarkerSize',8,'MarkerFaceColor',colors_plot(mm,:))
hold on
        legendnames{1,size(legendnames,2)+1} = fieldnames_beh{mm+1};
   end
end

for mm =1:max((candidate_frames))
   [ framesubset_cand,framesubset_tsne,~] = intersect(subset_of_points_to_plot_emreal,find(candidate_frames == mm));

   if numel(framesubset_cand)
figure(559)
plot3(ML_features.pose_score(framesubset_cand,1),ML_features.pose_score(framesubset_cand,2),...
   ML_features.pose_score(framesubset_cand,3),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))

% plot3(ML_features.ja_dyadic_spectrograms(framesubset_cand,1),ML_features.ja_dyadic_spectrograms(framesubset_cand,2),...
%     ML_features.ja_dyadic_spectrograms(framesubset_cand,3),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))

%plot3(ML_features.spectrogram_pcs_trunk_angle(framesubset_cand,1),ML_features.spectrogram_pcs_trunk_angle(framesubset_cand,2),...
 %   ML_features.spectrogram_pcs_trunk_angle(framesubset_cand,3),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))


%ML_features.ja_dyadic_spectrograms;

figure(572)
%plot(mappedX_dyn_angle(framesubset_tsne,1),mappedX_dyn_angle(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
%plot(mappedX_dyn_head(framesubset_tsne,1),mappedX_dyn_head(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
plot( zValues(framesubset_tsne,1), zValues(framesubset_tsne,2),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))


   end
end

 figure(559)
legend(legendnames)

 figure(572)
legend(legendnames)

make_ethogram(outputvector,fieldnames_beh)
        make_dotplot(    candidate_frames,fieldnames_beh,ML_features)
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
