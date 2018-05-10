function behavioral_error_output_struct = spec_beh_error(ML_features,cp,frames,labels,cp_num,fieldnames_beh,chosen_behavior)
    
    % create training and testing sets
    training_frames = frames(training(cp,cp_num));
    training_labels = labels(training(cp,cp_num));
    testing_frames = frames(test(cp,cp_num));
    
    % find number for desired behavior
    behavior_number = find(strcmp(fieldnames_beh,chosen_behavior)==1) - 1;
   
    % train model on training set, create prediction for testing set
    [treebagger_labels,~,~] = ...
    findsimilarframes_mlfeatures_multifeatures(ML_features,training_labels,training_frames,testing_frames);
    
    % initialize values
    num_annotated = 0;
    num_predicted = 0;
    num_correct_pred = 0;
    
    % calculate various desired outputs
    hand_labels = labels(test(cp,cp_num));
    for n = 1:length(testing_frames)
       if hand_labels(n) == behavior_number
           num_annotated = num_annotated + 1;
           if treebagger_labels(n) == behavior_number
               num_predicted = num_predicted + 1; 
               num_correct_pred = num_correct_pred + 1;
           end
       elseif treebagger_labels(n) == behavior_number
           num_predicted = num_predicted + 1;
       end
    end
    
    % save output to struct
    behavioral_error_output_struct.num_annotated = num_annotated;
    behavioral_error_output_struct.num_predicted = num_predicted;
    behavioral_error_output_struct.num_correct_pred = num_correct_pred;
    behavioral_error_output_struct.num_fp = num_predicted - num_correct_pred;
    behavioral_error_output_struct.num_fn = num_annotated - num_correct_pred;

end