%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATES CROSS-VALIDATED ERROR IN VARIOUS FORMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use this to save relevant workspace variables in a file called 'data.mat'
save('data2','frames','labels','ML_features','fieldnames_beh','observed_behaviors')

% rename annotated frames and corresponding labels
frames = annotated_frames_goodtracking;
labels = annotation_output_goodtracking;

% create partition
num_folds = 12;
cp = cvpartition(frames,'kfold',num_folds);

%% OVERALL CV ERROR

% run treebagger
cv_errors = zeros(1,num_folds);
confusmat = zeros(length(unique(labels)));
cfmat = zeros(length(fieldnames_beh)-1);
parfor cvindex = 1:num_folds
    [current_error,current_confusmat,current_cfmat] = cv_error(ML_features,cp,frames,labels,cvindex,fieldnames_beh);
    cv_errors(cvindex) = current_error;
    confusmat = confusmat + current_confusmat;
    cfmat = cfmat + current_cfmat;
end   

% average cv error
cv_error_avg = mean(cv_errors);

% convert cfmat to percentage
for i = 1:length(cfmat)
    cfmat(i,:) = cfmat(i,:) / sum(cfmat(i,:));
end

% invert matrix to make for more readable output to excel
cfmat = cfmat';

%% BEHAVIOR SPECIFIC ERROR (for one behavior only)

chosen_behavior = 'Walk';

num_annotated = 0;
num_predicted = 0;
num_correct_pred = 0;

% run RF on each training set, generate prediction for each corresponding
% testing set
parfor behaviorindex = 1:num_folds
    behavioral_output = spec_beh_error(ML_features,cp,frames,labels,behaviorindex,fieldnames_beh,chosen_behavior);
    num_annotated = num_annotated + behavioral_output.num_annotated;
    num_predicted = num_predicted + behavioral_output.num_predicted;
    num_correct_pred = num_correct_pred + behavioral_output.num_correct_pred;
end

% save output to struct
behavioral_error.num_annotated = num_annotated;
behavioral_error.num_predicted = num_predicted;
behavioral_error.num_correct_pred = num_correct_pred;
behavioral_error.num_fp = num_predicted - num_correct_pred;
behavioral_error.num_fn = num_annotated - num_correct_pred;

% misclassification rate for specified behavior 
mcrate = 1 - (behavioral_error.num_correct_pred / behavioral_error.num_annotated);

%%  COMPLETE BEHAVIORAL ERROR (outputs one struct with all behaviors)

cfmat = zeros(length(fieldnames_beh) - 1);
% one struct for each fold of partition
parfor index = 1:num_folds
    [beh_folds(index),current_cfmat] = spec_beh_error2(ML_features,cp,frames,labels,index,fieldnames_beh,fieldnames_beh);
    cfmat = cfmat + current_cfmat ; 
end

fields = {'num_annotated', 'num_predicted', 'num_correct_pred'};

% iterate through structs, creating the desired final struct 
for i = 2:numel(fieldnames_beh)   
    for j = 1:numel(fields)
        all_behaviors.(fieldnames_beh{i}).(fields{j}) = 0;        
        for k = 1:length(beh_folds)
            all_behaviors.(fieldnames_beh{i}).(fields{j}) = all_behaviors.(fieldnames_beh{i}).(fields{j}) + beh_folds(k).(fieldnames_beh{i}).(fields{j});
        end
    end     
    all_behaviors.(fieldnames_beh{i}).num_fp = all_behaviors.(fieldnames_beh{i}).num_predicted - all_behaviors.(fieldnames_beh{i}).num_correct_pred ;
    all_behaviors.(fieldnames_beh{i}).num_fn = all_behaviors.(fieldnames_beh{i}).num_annotated - all_behaviors.(fieldnames_beh{i}).num_correct_pred ;
    all_behaviors.(fieldnames_beh{i}).mc_rate = 1 - (all_behaviors.(fieldnames_beh{i}).num_correct_pred / all_behaviors.(fieldnames_beh{i}).num_annotated) ;
    all_behaviors.(fieldnames_beh{i}).fp_rate = all_behaviors.(fieldnames_beh{i}).num_fp / all_behaviors.(fieldnames_beh{i}).num_predicted ;
    all_behaviors.(fieldnames_beh{i}).fn_rate = all_behaviors.(fieldnames_beh{i}).num_fn / all_behaviors.(fieldnames_beh{i}).num_annotated ;
    all_behaviors.(fieldnames_beh{i}).prior_prob = all_behaviors.(fieldnames_beh{i}).num_annotated / length(frames);
    all_behaviors.(fieldnames_beh{i}).post_prob = all_behaviors.(fieldnames_beh{i}).num_predicted / length(frames);
    all_behaviors.(fieldnames_beh{i}).true_pos_rate = all_behaviors.(fieldnames_beh{i}).num_correct_pred / all_behaviors.(fieldnames_beh{i}).num_annotated ;
end
    
    % CONFUSION MATRIX 
    % convert cfmat to percentage
    for i = 1:length(cfmat)
        cfmat(i,:) = cfmat(i,:) / sum(cfmat(i,:));
    end
    % invert confusion matrix to make for more readable output to excel
    cfmat = cfmat';

    % TRUE POSITIVE VECTOR 
    behs = fieldnames(all_behaviors);
    tp_rate = [] ; 
    for i = 1:numel(behs)  
      if all_behaviors.(behs{i}).num_annotated ~= 0   
         tp_rate = [tp_rate, all_behaviors.(behs{i}).true_pos_rate]; 
      end
    end
    tp_rate = tp_rate' ; 
    
    % NUM_ANNOTATED VECTOR 
    behs = fieldnames(all_behaviors);
    num_annotated = [] ; 
    for i = 1:numel(behs)  
      if all_behaviors.(behs{i}).num_annotated ~= 0   
         num_annotated = [num_annotated, all_behaviors.(behs{i}).num_annotated]; 
      end
    end
    num_annotated = num_annotated' ;

% OUTPUT
% To display all of the fields in command window at once
behs = fieldnames(all_behaviors);
for i = 1:numel(behs)  
  if all_behaviors.(behs{i}).num_annotated ~= 0   
      fprintf(behs{i})
      fprintf('\n')
      display(all_behaviors.(behs{i})) 
  end
end

% outputs
all_behaviors ;
cfmat ;
tp_rate ;
num_annotated ; 

%% BLOCKED CV ERROR - (one fold)

% train on first 90%, test on last 10%  
cutoff = floor(0.95 * length(frames));

% create training and testing sets
training_frames = frames(1:cutoff);
training_labels = labels(1:cutoff);
testing_frames = frames((cutoff + 1):end);

% train model on training set, create prediction for testing set
[treebagger_labels,~,~] = ...
findsimilarframes_mlfeatures_multifeatures(ML_features,training_labels,training_frames,testing_frames);

% calculate cross validation error
num_correct = 0;
hand_labels = labels((cutoff + 1):end);
for n = 1:length(hand_labels)
   if hand_labels(n) == treebagger_labels(n)
       num_correct = num_correct + 1;
   end
end
blocked_cv_error = 1 - (num_correct / length(testing_frames));

%% BLOCKED CV ERROR - K-FOLD

num_folds = 3;

% need to drop out up to num_folds number of frames off of end because reshape must produce a valid matrix
num_to_drop = mod(length(frames),num_folds);
% create the blocked partition
frames_partition = reshape(frames(1:(end - num_to_drop)),[],num_folds);
labels_partition = reshape(labels(1:(end - num_to_drop)),[],num_folds);
% initialization
indices = 1:num_folds;
blocked_cv_errors = zeros(1,num_folds);
blocked_cfmat = zeros(length(fieldnames_beh) - 1);

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
    
    % confusion matrix for this fold
    current_cfmat = zeros(length(fieldnames_beh) - 1);
    for n = 1:(length(fieldnames_beh) - 1)
        for k = 1:(length(hand_labels))
            if hand_labels(k) == n
                num = treebagger_labels(k);
                current_cfmat(n,num) = current_cfmat(n,num) + 1;
            end
        end
    end 
    % add cfmat for this fold to overall cfmat
    blocked_cfmat = blocked_cfmat + current_cfmat;  
end

% compute average blocked cv error
blocked_cv_error_avg = mean(blocked_cv_errors);

% convert cfmat to percentage
for i = 1:length(blocked_cfmat)
    blocked_cfmat(i,:) = blocked_cfmat(i,:) / sum(blocked_cfmat(i,:));
end

% invert matrix to make for more readable output to excel
blocked_cfmat = blocked_cfmat';

%% BEHAVIORAL ERROR FOR BLOCKED PARTITION - (one behavior)

chosen_behavior = 'Walk';

% train on first 90%, test on last 10%  
cutoff = floor(0.9 * length(frames));

% create training and testing sets
training_frames = frames(1:cutoff);
training_labels = labels(1:cutoff);
testing_frames = frames((cutoff + 1):end);

% run treebagger
[treebagger_labels,~,~] = ...
findsimilarframes_mlfeatures_multifeatures(ML_features,training_labels,training_frames,testing_frames);

% initialize values
num_annotated = 0;
num_predicted = 0;
num_correct_pred = 0;

% find number for desired behavior
behavior_number = find(strcmp(fieldnames_beh,chosen_behavior)==1) - 1;
    
% calculate various desired outputs
hand_labels = labels((cutoff + 1):end);
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
blocked_behavioral_error.num_annotated = num_annotated;
blocked_behavioral_error.num_predicted = num_predicted;
blocked_behavioral_error.num_correct_pred = num_correct_pred;
blocked_behavioral_error.num_fp = num_predicted - num_correct_pred;
blocked_behavioral_error.num_fn = num_annotated - num_correct_pred;

%% BEHAVIORAL ERROR, CV ERROR, AND CONFUSMAT ON BLOCKED PARTITION - K-FOLD

    num_folds = 12;

    % need to drop out up to num_folds number of frames off of end because reshape must produce a valid matrix
    num_to_drop = mod(length(frames),num_folds);
    % create the blocked partition
    frames_partition = reshape(frames(1:(end - num_to_drop)),[],num_folds);
    labels_partition = reshape(labels(1:(end - num_to_drop)),[],num_folds);
    % initialization
    indices = 1:num_folds;
    blocked_cv_errors = zeros(1,num_folds);
    blocked_cfmat = zeros(length(fieldnames_beh) - 1);

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
        % create vector for testing set, corresponding labels
        testing_frames = frames_partition(:,i);  
        hand_labels = labels_partition(:,i);
        % run treebagger
%         [treebagger_labels,~] = ...
%         findsimilarframes_mlfeatures_multifeatures(ML_features,training_labels,training_frames,testing_frames);
        classnames = [] ; 
        cost_mat = [] ; 
        [treebagger_labels,~] = ...
        findsimilarframes_mlfeatures_multifeatures_2(ML_features,training_labels,training_frames,testing_frames,classnames,cost_mat);

        % BEHAVIORAL OUTPUT (save this fold's output to struct array)  
        beh_folds(i) = spec_beh_error3(fieldnames_beh,fieldnames_beh,hand_labels,treebagger_labels);

        % CV ERROR (for this fold)
        num_correct = 0;
        hand_labels = labels_partition(:,i);
        for n = 1:length(testing_frames)
           if hand_labels(n) == treebagger_labels(n)
               num_correct = num_correct + 1;
           end
        end
        % add cv error for this run to array
        blocked_cv_errors(i) = 1 - (num_correct / length(testing_frames));

        % CONFUSION MATRIX for this fold
        current_cfmat = zeros(length(fieldnames_beh) - 1);
        for n = 1:(length(fieldnames_beh) - 1)
            for k = 1:(length(hand_labels))
                if hand_labels(k) == n
                    num = treebagger_labels(k);
                    current_cfmat(n,num) = current_cfmat(n,num) + 1;
                end
            end
        end 
        % add cfmat for this fold to overall cfmat
        blocked_cfmat = blocked_cfmat + current_cfmat;  
    end

    % CV ERROR
    % compute average blocked cv error
    blocked_cv_error_avg = mean(blocked_cv_errors);

    % CONFUSION MATRIX
    % convert cfmat to percentage
    for i = 1:length(blocked_cfmat)
        blocked_cfmat(i,:) = blocked_cfmat(i,:) / sum(blocked_cfmat(i,:));
    end
    % invert confusion matrix to make for more readable output to excel
    blocked_cfmat = blocked_cfmat';

    % BEHAVIORAL OUTPUT
    fields = {'num_annotated', 'num_predicted', 'num_correct_pred'};
    % iterate through structs, creating the desired final struct 
    for i = 2:numel(fieldnames_beh)   
        for j = 1:numel(fields)
            all_behaviors.(fieldnames_beh{i}).(fields{j}) = 0;        
            for k = 1:length(beh_folds)
                all_behaviors.(fieldnames_beh{i}).(fields{j}) = all_behaviors.(fieldnames_beh{i}).(fields{j}) + beh_folds(k).(fieldnames_beh{i}).(fields{j});
            end
        end     
        all_behaviors.(fieldnames_beh{i}).num_fp = all_behaviors.(fieldnames_beh{i}).num_predicted - all_behaviors.(fieldnames_beh{i}).num_correct_pred ;
        all_behaviors.(fieldnames_beh{i}).num_fn = all_behaviors.(fieldnames_beh{i}).num_annotated - all_behaviors.(fieldnames_beh{i}).num_correct_pred ;
        all_behaviors.(fieldnames_beh{i}).mc_rate = 1 - (all_behaviors.(fieldnames_beh{i}).num_correct_pred / all_behaviors.(fieldnames_beh{i}).num_annotated) ;
        all_behaviors.(fieldnames_beh{i}).fp_rate = all_behaviors.(fieldnames_beh{i}).num_fp / all_behaviors.(fieldnames_beh{i}).num_predicted ;
        all_behaviors.(fieldnames_beh{i}).fn_rate = all_behaviors.(fieldnames_beh{i}).num_fn / all_behaviors.(fieldnames_beh{i}).num_annotated ;
        all_behaviors.(fieldnames_beh{i}).prior_prob = all_behaviors.(fieldnames_beh{i}).num_annotated / (length(frames) - num_to_drop);
        all_behaviors.(fieldnames_beh{i}).post_prob = all_behaviors.(fieldnames_beh{i}).num_predicted / (length(frames) - num_to_drop);
        all_behaviors.(fieldnames_beh{i}).true_pos_rate = all_behaviors.(fieldnames_beh{i}).num_correct_pred / all_behaviors.(fieldnames_beh{i}).num_annotated ;
    end
    
    % TRUE POSITIVE VECTOR 
    behs = fieldnames(all_behaviors);
    tp_rate = [] ; 
    for i = 1:numel(behs)  
      if all_behaviors.(behs{i}).num_annotated ~= 0   
         tp_rate = [tp_rate, all_behaviors.(behs{i}).true_pos_rate]; 
      end
    end
    tp_rate = tp_rate' ; 
    
    % NUM_ANNOTATED VECTOR 
    behs = fieldnames(all_behaviors);
    num_annotated = [] ; 
    for i = 1:numel(behs)  
      if all_behaviors.(behs{i}).num_annotated ~= 0   
         num_annotated = [num_annotated, all_behaviors.(behs{i}).num_annotated]; 
      end
    end
    num_annotated = num_annotated' ;
    
    

% OUTPUT
% To display all of the behavioral fields in command window at once
behs = fieldnames(all_behaviors);
for i = 1:numel(behs)  
  if all_behaviors.(behs{i}).num_annotated ~= 0   
      fprintf(behs{i})
      fprintf('\n')
      display(all_behaviors.(behs{i})) 
  end
end
   
% outputs
all_behaviors ;
blocked_cfmat ;
blocked_cv_errors ;
blocked_cv_error_avg ;
tp_rate ;
num_annotated ; 

%% Instanced Behavior

num_folds = 8;

    

    % need to drop out up to num_folds number of frames off of end because reshape must produce a valid matrix
    num_to_drop = mod(length(frames),num_folds);
    % create the blocked partition
    frames_partition = reshape(frames(1:(end - num_to_drop)),[],num_folds);
    labels_partition = reshape(labels(1:(end - num_to_drop)),[],num_folds);
    % initialization
    indices = 1:num_folds;
    blocked_cv_errors = zeros(1,num_folds);
    blocked_cfmat = zeros(length(fieldnames_beh));

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
        
        training_labels = training_labels(:) ; 
        training_frames = training_frames(:) ; 
        
        testing_frames = frames_partition(:,i); 
        hand_labels = labels_partition(:,i);   
        
        % instance the training labels and features
        [instanced_labels,instanced_features,Mdl_test] = instance_data(ML_features,training_labels,training_frames);
        
        [instanced_labels_test,instanced_features_test,~] = instance_data(ML_features,hand_labels,testing_frames);
        
       
        % run prediction for heldout FRAMES (not instances)
        [testval,~] = Mdl_test.predict(instanced_features_test);
        treebagger_labels = cellfun(@str2num,testval);

        % BEHAVIORAL OUTPUT (save this fold's output to struct array)  
%         beh_folds(i) = spec_beh_error3(fieldnames_beh,observed_behaviors,hand_labels,treebagger_labels);
        [beh_folds(i),~] = beh_analysis_3(instanced_labels_test, treebagger_labels, fieldnames_beh)

        % CV ERROR (for this fold)
        num_correct = 0;
        
        for n = 1:length(instanced_labels_test)
           if instanced_labels_test(n) == treebagger_labels(n)
               num_correct = num_correct + 1;
           end
        end
        % add cv error for this run to array
        blocked_cv_errors(i) = 1 - (num_correct / length(testing_frames));

        % CONFUSION MATRIX for this fold
        current_cfmat = zeros(length(fieldnames_beh));
        for n = 1:(length(fieldnames_beh))
            for k = 1:(length(instanced_labels_test))
                if instanced_labels_test(k) == n
                    num = treebagger_labels(k);
                    current_cfmat(n,num) = current_cfmat(n,num) + 1;
                end
            end
        end 
        % add cfmat for this fold to overall cfmat
        blocked_cfmat = blocked_cfmat + current_cfmat;  
    end

    % CV ERROR
    % compute average blocked cv error
    blocked_cv_error_avg = mean(blocked_cv_errors);

    % CONFUSION MATRIX
    % convert cfmat to percentage
    for i = 1:length(blocked_cfmat)
        blocked_cfmat(i,:) = blocked_cfmat(i,:) / sum(blocked_cfmat(i,:));
    end
    % invert confusion matrix to make for more readable output to excel
    blocked_cfmat = blocked_cfmat';

    % BEHAVIORAL OUTPUT
    fields = {'num_annotated', 'num_predicted', 'num_correct_pred'};
    % iterate through structs, creating the desired final struct 
    for i = 1:length(fieldnames_beh)   
        for j = 1:numel(fields)
            all_behaviors.(fieldnames_beh{i}).(fields{j}) = 0;        
            for k = 1:length(beh_folds)
                all_behaviors.(fieldnames_beh{i}).(fields{j}) = all_behaviors.(fieldnames_beh{i}).(fields{j}) + beh_folds(k).(fieldnames_beh{i}).(fields{j});
            end
        end     
        all_behaviors.(fieldnames_beh{i}).num_fp = all_behaviors.(fieldnames_beh{i}).num_predicted - all_behaviors.(fieldnames_beh{i}).num_correct_pred ;
        all_behaviors.(fieldnames_beh{i}).num_fn = all_behaviors.(fieldnames_beh{i}).num_annotated - all_behaviors.(fieldnames_beh{i}).num_correct_pred ;
        all_behaviors.(fieldnames_beh{i}).mc_rate = 1 - (all_behaviors.(fieldnames_beh{i}).num_correct_pred / all_behaviors.(fieldnames_beh{i}).num_annotated) ;
        all_behaviors.(fieldnames_beh{i}).fp_rate = all_behaviors.(fieldnames_beh{i}).num_fp / all_behaviors.(fieldnames_beh{i}).num_predicted ;
        all_behaviors.(fieldnames_beh{i}).fn_rate = all_behaviors.(fieldnames_beh{i}).num_fn / all_behaviors.(fieldnames_beh{i}).num_annotated ;
        all_behaviors.(fieldnames_beh{i}).prior_prob = all_behaviors.(fieldnames_beh{i}).num_annotated / (length(frames) - num_to_drop);
        all_behaviors.(fieldnames_beh{i}).post_prob = all_behaviors.(fieldnames_beh{i}).num_predicted / (length(frames) - num_to_drop);
        all_behaviors.(fieldnames_beh{i}).true_pos_rate = all_behaviors.(fieldnames_beh{i}).num_correct_pred / all_behaviors.(fieldnames_beh{i}).num_annotated ;
    end
    
    % TRUE POSITIVE VECTOR 
    behs = fieldnames(all_behaviors);
    tp_rate = [] ; 
    for i = 1:numel(behs)  
      if all_behaviors.(behs{i}).num_annotated ~= 0   
         tp_rate = [tp_rate, all_behaviors.(behs{i}).true_pos_rate]; 
      end
    end
    tp_rate = tp_rate' ; 
    
    % NUM_ANNOTATED VECTOR 
    behs = fieldnames(all_behaviors);
    num_annotated = [] ; 
    for i = 1:numel(behs)  
      if all_behaviors.(behs{i}).num_annotated ~= 0   
         num_annotated = [num_annotated, all_behaviors.(behs{i}).num_annotated]; 
      end
    end
    num_annotated = num_annotated' ;
      
    

% OUTPUT
% To display all of the behavioral fields in command window at once
behs = fieldnames(all_behaviors);
for i = 1:numel(behs)  
  if all_behaviors.(behs{i}).num_annotated ~= 0   
      fprintf(behs{i})
      fprintf('\n')
      display(all_behaviors.(behs{i})) 
  end
end
   
% outputs
all_behaviors ;
blocked_cfmat ;
blocked_cv_errors ;
blocked_cv_error_avg ;
tp_rate ;
num_annotated ; 

% create_heatmap( CFMAT_PERCENT, CFMAT_RAWNUM, FIELDNAMES, SWITCH FLAG, FIGURE NUMBER)
create_heatmap(blocked_cfmat, [], fieldnames_beh, 1,1 ) ;
title('Blocked CV Error 8fold Vicon8 Caff')




changes = []; 
for i = 1:(length(labels) - 1)
    if labels(i) ~= labels(i + 1) 
        changes = [changes ; labels(i + 1)] ; 
    end
end





