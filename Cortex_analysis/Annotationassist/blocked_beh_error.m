function behavioral_error_output_struct = blocked_beh_error(ML_features,frames,labels,fieldnames_beh,chosen_behavior)
    
    % need to drop out up to num_folds number of frames off of end 
    % because reshape must produce a valid matrix
    num_to_drop = mod(length(frames),num_folds);
    frames = frames(1:(end - num_to_drop));
    labels = labels(1:(end - num_to_drop));
    
    % create the blocked partition
    frames_partition = reshape(frames,[],num_folds);
    labels_partition = reshape(labels,[],num_folds);

    % initialization
    indices = 1:num_folds;
    blocked_cv_errors = zeros(1,num_folds);

    for i = 1:num_folds
    
        % determine which part to hold out
        training_items = indices(indices~=i);

        % initialization
        training_labels = [0];
        training_frames = [0];

        % create vector for training set
        for j = 1:(num_folds - 1)
            j2 = training_items(j);
            training_labels = [training_labels,labels_partition(:,j2)'];
            training_frames = [training_frames,frames_partition(:,j2)'];              
        end

        % create vector for testing set
        testing_frames = frames_partition(:,i);

        % run treebagger
        [treebagger_labels,~,~] = ...
        findsimilarframes_mlfeatures_multifeatures(ML_features,training_labels,training_frames,testing_frames);
    end 
    
    % initialize values
    num_annotated = 0;
    num_predicted = 0;
    num_correct_pred = 0;
    
    % calculate various desired outputs
    hand_labels = labels_partition(:,i);
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