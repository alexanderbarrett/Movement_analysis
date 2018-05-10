%% set this to the server address
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';
%mocapmasterdirectory = '\\10.242.62.28\olveczky_lab\Jesse\Data\Motionanalysis_captures\';

%% PARMAETERS
%for the ML features
ratname = 'Vicon8';
conditionnumbers = 6;
annotation_number = 2;

% JDM25
ratname = 'JDM25';
conditionnumbers = 5; %1 is pre lesion, 10 is post lesion
 
% JDM33
% ratname = 'JDM33';
% conditionnumber = [1,5,10]; %1 is pre lesion, 10 is post lesion

 
% annotation_number = 4; % look in the file below to check which file to load in
 
[fullposture_annotation_struct,fulloutput_annotation_struct,subset_of_frames_annotated] = load_handannotation_files(annotation_number,mocapmasterdirectory);
 
%% get files and the machine learning features
descriptor_struct = get_mocap_files_table(conditionnumbers(1),ratname,mocapmasterdirectory);
[~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
%[mocapstruct_all] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,0,[],mocapvideodirectory,1);
[~,mocapstruct_matobjs] = preprocess_mocap_data_3(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,0,[],mocapvideodirectory,0,1);
mocapstruct_all = load(mocapstruct_matobjs{1}.Properties.Source);
for ll = 2:numel(mocapstruct_matobjs)
mocapstruct_all =mergemocapstructs(mocapstruct_all,load(mocapstruct_matobjs{ll}.Properties.Source));
end
 
annotation_string = 'annotation_restricted_ii';
%%%%%%%% NOTE - took out overwriteML as third argument here %%%%%%%%%%%%%%%
ML_matobjs{1} = batch_compute_mlfeatures(mocapfilearray,ratname,0);
 
%% load in the ML features
[ML_features] = load_selectML(1,ML_matobjs,annotation_string);
   % load_tsne_mocap(1,ML_matobjs,mocapstruct_matobjs,annotation_string,1,'fast');
 
%ML_features = get_super vised_features(mocapstruct_all,mocapstruct_all.modular_cluster_properties.clustering_inds_agg{2},2,ratname,descriptor_struct.Nametag,0,0);
 
%% visualize a chosen behavior (optional)
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
 
 
 
%% RUN THE ML ALGORITHM TO FIND SIMILAR FRAMES
fprintf('searching for similar frames \n')
frames_to_annotate = 1:max(subset_of_frames_annotated);
[candidate_frames,score,Mdl_test] = findsimilarframes_mlfeatures_multifeatures_2(ML_features,annotation_output_goodtracking,annotated_frames_goodtracking,frames_to_annotate);

oob_err = oobError(Mdl_test,'mode','ensemble'); 
 
 
%% VISUALIZE AN EXAMPLE BEHAVIOR

chosenbehavior = 'WetDogShake';
animate_markers_aligned_fullmovie(mocapstruct_all,fulloutput_annotation_struct.(chosenbehavior)(1:10:end)')

 
candidate_frames_vec = (find(candidate_frames==find(strcmp(fieldnames_beh ,chosenbehavior)==1)));
candidate_frame_predicted = setxor(ML_features.framelist_true(candidate_frames_vec),...
    intersect(ML_features.framelist_true,fulloutput_annotation_struct.(chosenbehavior)));
 
%visualize the ground truth and predicted data
animate_markers_aligned_fullmovie(mocapstruct_all,(candidate_frame_predicted(1:10:end))')
 
 
 
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