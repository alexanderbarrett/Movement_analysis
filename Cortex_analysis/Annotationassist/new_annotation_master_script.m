
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';

%% PARMAETERS
%for the ML features
ratname = 'Vicon8';
conditionname = 'Vicon8_caff';
conditionnumber = 6;

% for the specific file
chosenbehavior = 'WetDogShake';
subset_of_frames_to_annotate = 1:2*10^6;
gapfill_number = 20;

%% can load in multiple filenames
annotation_mocapname = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Preprocessed\nolj_Recording_day7_caffeine1_nolj.mat';
annotation_filename = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Preprocessed\nolj_Recording_day7_caffeine1_nolj_handannotation_JDM.mat';
annotation_filenumber = 1;
outputstruct_annotation = load(annotation_filename);
subset_of_frames_annotated = (annotation_filenumber-1)*540000+1:100000; %depends on the file(s) loaded


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
[fulloutput_annotation_struct,indivbouts_caff_struct] = fillannotationgaps_struct(outputstruct_annotation.output.GlobalBehavior,gapfill_number);

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
[candidate_frames,score,Mdl_test] = findsimilarframes_mlfeatures_multifeatures(ML_features,annotation_output_goodtracking,annotated_frames_goodtracking,subset_of_frames_to_annotate);
%get only the types of behavior that are found
unique_predictions = unique(candidate_frames(~isnan(candidate_frames)));
%these are the fieldnames of these behaviors
fieldnames_beh(unique_predictions+1)

%% VISUALIZE AN EXAMPLE BEHAVIOR
candidate_frames_vec = (find(candidate_frames==find(strcmp(fieldnames_beh ,'WetDogShake')==1)-1));
predicted_annotatedframes= candidate_frames_vec(candidate_frames_vec<80000);
predicted_unannotatedframes= candidate_frames_vec(candidate_frames_vec>80000);

%visualize the ground truth and predicted data
animate_markers_aligned_fullmovie(mocapstruct_all,frames_with_goodtracking(predicted_unannotatedframes(1:10:end))')
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


