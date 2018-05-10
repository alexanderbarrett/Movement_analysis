%% PARMAETERS

% initialize
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';
ratnames_master = ["Vicon8";"Vicon8";"JDM25";"JDM33";"JDM25";"JDM25";"JDM33";"JDM33"];
condition_numbers_vec_master = [6;1;5;1;1;1;1;1];
annotation_numbers_master = [2;3;4;5;7;8;9;10];

% select rats to load
    % 1 is Vicon8 caff, 2 is Vicon8 prelesion, 3 is JDM25 caff, 4 is JDM33 caff, 
    % 5 is JDM25 prelesion, 6 is JDM25 postlesion, 9 is JDM33 postlesion, 10 is JDM33 postlesion2
indices = [1; 3];
ratnames = ratnames_master(indices) ;
condition_numbers_vec = condition_numbers_vec_master(indices) ;
annotation_numbers = annotation_numbers_master(indices) ;

% load rats
[tot_ML_features, tot_fieldnames_beh, tot_annotation_output_goodtracking, tot_annotated_frames_goodtracking] ...
    = load_rat_data(ratnames, condition_numbers_vec, annotation_numbers, mocapmasterdirectory, struct([]), {}, [], [], []) ; 

% save a copy / rename
frames = tot_annotated_frames_goodtracking ;
labels = tot_annotation_output_goodtracking ; 
fieldnames_beh = tot_fieldnames_beh ; 
features = tot_ML_features ; 

% equivalence classes
[labels,fieldnames_beh] = get_equivalence_classes(labels,fieldnames_beh) ;

% remove unwanted behaviors
% bad_behs = ["Swap","Transition","HeadSwap","FSpineSwap","FspineSwap"];
bad_behs = ["Transition"] ; 
for i = 1:length(bad_behs)
    current_beh = char(bad_behs(i)) ;
    bad_num = find(strcmp(fieldnames_beh,current_beh)==1);
    indices = labels ~= bad_num;
    frames = frames(indices);
    labels = labels(indices);
end

%% RUN TREEBAGGER ON ITSELF (cross validation)

[candidate_frames,score,Mdl_test,annotated_instances] = findsimilarframes_mlfeatures_multifeatures(features,labels,frames,[]);
[all_behaviors,cfmat_percent,cfmat_rawnum,tp_rate] = beh_analysis_3(annotated_instances,candidate_frames,fieldnames_beh);
% create heatmap (1 is percent, 2 is raw nums)
create_heatmap(cfmat_percent, cfmat_rawnum, fieldnames_beh, 1, 1) ;
figure(1) 
title('Cross Validation (30% heldout) : JDM25 caff and Vicon8 caff');


%% TEST NEW RATS

% load desired model
trained_treebagger = load('Vicon8caff_Mdltest_annotation_ii.mat') ;
trained_Mdl_test = trained_treebagger.Mdl_test ; 
Reference_fieldnames = load('Vicon8caff_fixedtotfieldnamesbeh_annotation_ii.mat') ;
good_fieldnames = Reference_fieldnames.tot_fieldnames_beh ;

% match new rat's labels to model test
[good_fieldnames,labels] = adjust_label_nums(fieldnames_beh,good_fieldnames,labels) ; 

% test new rat in model
[candidate_frames,annotated_instances] = test_new_rat(features,labels,frames,trained_Mdl_test);
[all_behaviors,cfmat_percent,cfmat_rawnum,tp_rate] = beh_analysis_3(annotated_instances,candidate_frames,good_fieldnames);

% create heatmap (3 is percent, 4 is raw nums)
create_heatmap(cfmat_percent, cfmat_rawnum, good_fieldnames,4) ;
figure(2) 
title('testing');

%% display

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
percent = corr / ann ; 



