function cv_error_output = blocked_cv_error(ML_features,frames,labels,num_folds)
    
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

    parfor i = 1:num_folds
    
        % determine which part to hold out
        training_items = indices(indices~=i);

        % initialization
        training_labels = [];
        training_frames = [];

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

        % calculate cross validation error
        num_correct = 0;
        hand_labels = labels_partition(:,i);
        for n = 1:length(testing_frames)
           if hand_labels(n) == treebagger_labels(n)
               num_correct = num_correct + 1;
           end
        end
        blocked_cv_error = 1 - (num_correct / length(testing_frames));

        % add cv error for this run to array
        blocked_cv_errors(i) = blocked_cv_error;
    
    end
    
    % return error array
    cv_error_output = blocked_cv_errors;

end