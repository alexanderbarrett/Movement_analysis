
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';

%% PARMAETERS
%for the ML features
ratname = 'Vicon8';
conditionname = 'Vicon8_caff';
conditionnumber = 6;

% for the specific file
subset_of_frames_to_annotate = 1:2*10^6;
gapfill_number = 20;

%% can load in multiple filenames
annotation_mocapname = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Preprocessed\nolj_Recording_day7_caffeine1_nolj.mat';
annotation_filename = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Preprocessed\nolj_Recording_day7_caffeine1_nolj_handannotation_JDM.mat';
annotation_filenumber = 1;

annotation_folder = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Preprocessed\';
annotation_filenames = {'nolj_Recording_day7_caffeine1_nolj_handannotation_AB_1to48000.mat',...
    'nolj_Recording_day7_caffeine2_nolj_handannotation_AB_1to23800.mat',...
    'nolj_Recording_day7_caffeine2_nolj_handannotation_AB_23800to53323.mat',...
    'nolj_Recording_day7_caffeine2_nolj_handannotation_AB_53323to64603.mat'};
annotation_filenumber = {1,2,2,2};

annot_cell = cell(1,numel(annotation_filenames));
subset_of_frames_annotated = [];
agg_struct = struct([]);
for jj = 1:4
output1 = load(strcat(annotation_folder,annotation_filenames{jj}));
for kk = fieldnames(output1.output.GlobalBehavior)'
    output1.output.GlobalBehavior.(kk{1}) = output1.output.GlobalBehavior.(kk{1})+(annotation_filenumber{jj}-1)*540000;
end
annot_cell{jj} = output1;

  C = struct2cell(output1.output.GlobalBehavior);
minval = min([C{:}]);
maxval = max([C{:}]);

subset_of_frames_annotated = cat(2,subset_of_frames_annotated ,(annotation_filenumber{jj}-1)*540000+minval:maxval); %depends on the file(s) loaded
agg_struct = mergeStructs_JDM(agg_struct,output1.output.GlobalBehavior);

end
subset_of_frames_annotated = unique(subset_of_frames_annotated);


outputstruct_annotation = load(annotation_filename);


%% get files and the machine learning features
descriptor_struct = get_mocap_files_table(conditionnumber,ratname);
 [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
 [mocapstruct_all] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,0,[],mocapvideodirectory,0);
ML_features = get_supervised_features(mocapstruct_all,mocapstruct_all.modular_cluster_properties.clustering_inds_agg{2},2);
%the frames are the intersection of those that are returned from clipping
%and in the cluster of markers
frames_with_goodtracking = intersect(mocapstruct_all.modular_cluster_properties.clustering_inds_agg{2},mocapstruct_all.modular_cluster_properties.clipped_index{2});




%% LOAD IN THE ANNOTATION
% fill the gaps in the structure
% in the future can write code to accomodate multiple annotation_files 
[fulloutput_annotation_struct,indivbouts_annotation_struct] = fillannotationgaps_struct(agg_struct ,gapfill_number);

%visualize a chosen behavior (optional)
%animate_markers_aligned_fullmovie(mocapstruct_caff,fulloutput_caff_struct.(chosenbehavior)(1:10:end)')

% Convert the struct to a categorical output for the random forest
[outputvector,fieldnames_beh,observed_behaviors] = generate_categorical_output(fulloutput_annotation_struct,numel(subset_of_frames_annotated));
fieldnames_beh = cat(1,'null',fieldnames_beh); %there is a 0 for unannotated in the categorical vector

% get the outputvector
annotation_output_goodtracking = outputvector(intersect(frames_with_goodtracking,subset_of_frames_annotated));
% get the subset of overall frames with good tracking
[~,annotated_frames_goodtracking] = intersect(frames_with_goodtracking,subset_of_frames_annotated);

%% RUN THE ML ALGORITHM TO FIND SIMILAR FRAMES
fprintf('searching for similar frames \n')
[candidate_frames,score,Mdl_test] = findsimilarframes_mlfeatures_multifeatures(ML_features,annotation_output_goodtracking,annotated_frames_goodtracking,subset_of_frames_to_annotate);
%get only the types of behavior that are found
unique_predictions = unique(candidate_frames(~isnan(candidate_frames)));
%these are the fieldnames of these behaviors
fieldnames_beh(unique_predictions+1)

%% VISUALIZE AN EXAMPLE BEHAVIOR
chosenbehavior = 'FaceWipe';
animate_markers_aligned_fullmovie(mocapstruct_all,fulloutput_annotation_struct.(chosenbehavior)(1:10:end)')

candidate_frames_vec = (find(candidate_frames==find(strcmp(fieldnames_beh ,chosenbehavior)==1)-1));
candidate_frame_predicted = setxor(frames_with_goodtracking(candidate_frames_vec),...
    intersect(frames_with_goodtracking,fulloutput_annotation_struct.(chosenbehavior)));
%predicted_annotatedframes= candidate_frames_vec(candidate_frames_vec<80000);
%predicted_unannotatedframes= candidate_frames_vec(candidate_frames_vec>80000);

%visualize the ground truth and predicted data
animate_markers_aligned_fullmovie(mocapstruct_all,(candidate_frame_predicted(1:10:end))')
animate_markers_aligned_fullmovie(mocapstruct_all,frames_with_goodtracking(predicted_annotatedframes(1:10:end))')

%% add other functions to compare real and found

%visualize the behavioral segregation
%% plot on the tsne
%Need to refactor this
colors_plot = hsv(80);
colors_plot = colorcube(80)./2+lines(80)./2;
legendnames = cell(1,0);
subset_of_points_to_plot = 1:100:100*10000;

for mm =1:max(candidate_frames)
   [ framesubset_cand,framesubset_tsne,~] = intersect(subset_of_points_to_plot,find(candidate_frames == mm));
   if numel(framesubset_cand)
        legendnames{1,size(legendnames,2)+1} = fieldnames_beh{mm+1};

figure(559)
plot3(ML_features.pose_score(framesubset_cand,1),ML_features.pose_score(framesubset_cand,2),...
    ML_features.pose_score(framesubset_cand,3),'o','MarkerEdgeColor','none','MarkerSize',2,'MarkerFaceColor',colors_plot(mm,:))
hold on
   end
end

 figure(559)
legend(legendnames)

make_ethogram(outputvector,fieldnames_beh)
        make_dotplot(    candidate_frames,fieldnames_beh,ML_features)
plot_behavior_examples(mocapstruct_all,fulloutput_annotation_struct.GeneralGroom,indivbouts_annotation_struct.GeneralGroom)


