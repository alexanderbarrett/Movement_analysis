function [all_behaviors,current_cfmat] = spec_beh_error2(ML_features,cp,frames,labels,cp_num,fieldnames_beh,observed_behaviors)
    
    % create training and testing sets
    training_frames = frames(training(cp,cp_num));
    training_labels = labels(training(cp,cp_num));
    testing_frames = frames(test(cp,cp_num));
    hand_labels = labels(test(cp,cp_num));
       
    % call treebagger
    [treebagger_labels,~,~] = ...
    findsimilarframes_mlfeatures_multifeatures(ML_features,training_labels,training_frames,testing_frames);
    
    for index = 2:length(observed_behaviors)
        
        % get current behavior name
        chosen_behavior = observed_behaviors(index);
        
        % find number for desired behavior
        behavior_number = find(strcmp(fieldnames_beh,chosen_behavior)==1) - 1;
        
        % initialize values
        num_annotated = 0;
        num_predicted = 0;
        num_correct_pred = 0;
    
        % calculate various desired outputs  
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
        chosen_behavior_error_struct.num_annotated = num_annotated;
        chosen_behavior_error_struct.num_predicted = num_predicted;
        chosen_behavior_error_struct.num_correct_pred = num_correct_pred;
        
        all_behaviors.(observed_behaviors{index}) = chosen_behavior_error_struct;
        
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
        
    end
 
end