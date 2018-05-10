function [all_behaviors,current_cfmat_pct,cfmat_num,tp_rate] = beh_analysis(hand_labels, treebagger_labels, fieldnames_beh)
%     for index = 2:length(observed_behaviors)
%         
%         % get current behavior name
%         chosen_behavior = observed_behaviors(index);
%         
%         % find number for desired behavior
%         behavior_number = find(strcmp(fieldnames_beh,chosen_behavior)==1) - 1;
%         
%         % initialize values
%         num_annotated = 0;
%         num_predicted = 0;
%         num_correct_pred = 0;
%     
%         % calculate various desired outputs  
%         for n = 1:length(hand_labels)
%             if hand_labels(n) == behavior_number
%                 num_annotated = num_annotated + 1;
%                 if treebagger_labels(n) == behavior_number
%                     num_predicted = num_predicted + 1; 
%                     num_correct_pred = num_correct_pred + 1;
%                 end
%             elseif treebagger_labels(n) == behavior_number
%                 num_predicted = num_predicted + 1;
%             end
%         end
%         
%         i = index ; 
%         all_behaviors.(observed_behaviors{i}).num_annotated = num_annotated;
%         all_behaviors.(observed_behaviors{i}).num_predicted = num_predicted;
%         all_behaviors.(observed_behaviors{i}).num_correct_pred = num_correct_pred;     
%         all_behaviors.(observed_behaviors{i}).num_fp = all_behaviors.(observed_behaviors{i}).num_predicted - all_behaviors.(observed_behaviors{i}).num_correct_pred ;
%         all_behaviors.(observed_behaviors{i}).num_fn = all_behaviors.(observed_behaviors{i}).num_annotated - all_behaviors.(observed_behaviors{i}).num_correct_pred ;
%         all_behaviors.(observed_behaviors{i}).mc_rate = 1 - (all_behaviors.(observed_behaviors{i}).num_correct_pred / all_behaviors.(observed_behaviors{i}).num_annotated) ;
%         all_behaviors.(observed_behaviors{i}).fp_rate = all_behaviors.(observed_behaviors{i}).num_fp / all_behaviors.(observed_behaviors{i}).num_predicted ;
%         all_behaviors.(observed_behaviors{i}).fn_rate = all_behaviors.(observed_behaviors{i}).num_fn / all_behaviors.(observed_behaviors{i}).num_annotated ;
%         all_behaviors.(observed_behaviors{i}).prior_prob = all_behaviors.(observed_behaviors{i}).num_annotated / length(hand_labels);
%         all_behaviors.(observed_behaviors{i}).post_prob = all_behaviors.(observed_behaviors{i}).num_predicted / length(hand_labels);
%         all_behaviors.(observed_behaviors{i}).true_pos_rate = all_behaviors.(observed_behaviors{i}).num_correct_pred / all_behaviors.(observed_behaviors{i}).num_annotated ;
%           
%         % CONFUSION MATRIX 
%         current_cfmat = zeros(length(fieldnames_beh) - 1);
%         for n = 1:(length(fieldnames_beh) - 1)
%             for k = 1:(length(hand_labels))
%                 if hand_labels(k) == n
%                     num = treebagger_labels(k);
%                     current_cfmat(n,num) = current_cfmat(n,num) + 1;
%                 end
%             end
%         end
%         
%         % TRUE POSITIVE VECTOR 
%         behs = fieldnames(all_behaviors);
%         tp_rate = [] ; 
%         for i = 1:numel(behs)  
%           if all_behaviors.(behs{i}).num_annotated ~= 0   
%              tp_rate = [tp_rate, all_behaviors.(behs{i}).true_pos_rate]; 
%           end
%         end
%         tp_rate = tp_rate' ;
%         
%     end

    for index = 2:length(fieldnames_beh)
               
        behavior_number = index - 1 ; 
        
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
        
        i = index ; 
        all_behaviors.(fieldnames_beh{i}).num_annotated = num_annotated;
        all_behaviors.(fieldnames_beh{i}).num_predicted = num_predicted;
        all_behaviors.(fieldnames_beh{i}).num_correct_pred = num_correct_pred;     
        all_behaviors.(fieldnames_beh{i}).num_fp = num_predicted - num_correct_pred ;
        all_behaviors.(fieldnames_beh{i}).num_fn = num_annotated - num_correct_pred ;
        all_behaviors.(fieldnames_beh{i}).mc_rate = num_correct_pred / num_annotated ;
        all_behaviors.(fieldnames_beh{i}).fp_rate = all_behaviors.(fieldnames_beh{i}).num_fp / num_predicted ;
        all_behaviors.(fieldnames_beh{i}).fn_rate = all_behaviors.(fieldnames_beh{i}).num_fn / num_annotated ;
        all_behaviors.(fieldnames_beh{i}).prior_prob = num_annotated / length(hand_labels);
        all_behaviors.(fieldnames_beh{i}).post_prob = num_predicted / length(hand_labels);
        all_behaviors.(fieldnames_beh{i}).true_pos_rate = num_correct_pred / num_annotated ;
          
        
        
    end

        % CONFUSION MATRIX 
        current_cfmat = zeros(length(fieldnames_beh) - 1);
        for n = 1:(length(fieldnames_beh) - 1)
            for k = 1:(length(hand_labels))
                if hand_labels(k) == n
                    num = treebagger_labels(k) ;
                    current_cfmat(n,num) = current_cfmat(n,num) + 1;
                end
            end
        end
        cfmat_num = current_cfmat' ; 
        for i = 1:length(current_cfmat)
            current_cfmat(i,:) = current_cfmat(i,:) / sum(current_cfmat(i,:));
        end
        % invert confusion matrix to make for more readable output to excel
        current_cfmat_pct = current_cfmat';
        
        
        

        % TRUE POSITIVE VECTOR 
        behs = fieldnames(all_behaviors);
        tp_rate = [] ; 
        for i = 1:numel(behs)  
          if all_behaviors.(behs{i}).num_annotated ~= 0   
             tp_rate = [tp_rate, all_behaviors.(behs{i}).true_pos_rate]; 
          end
        end
        tp_rate = tp_rate' ;



end