function all_behaviors = spec_beh_error3(fieldnames_beh,observed_behaviors,hand_labels,treebagger_labels)
    
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
        for n = 1:length(hand_labels)
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
        
    end
 
end