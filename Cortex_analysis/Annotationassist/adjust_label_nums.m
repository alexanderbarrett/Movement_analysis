function [good_fieldnames_beh, tot_annotation_output_goodtracking] ...
    = adjust_label_nums(tot_fieldnames_beh, good_fieldnames_beh, tot_annotation_output_goodtracking)

original_labels_newfile = tot_annotation_output_goodtracking;
% change this rat's numbering system to match that of trained_treebagger
for j = 1:length(tot_fieldnames_beh)
    current_beh = char(tot_fieldnames_beh(j)) ; 
    % if current behavior has already been observed
    if sum(strcmp(current_beh,good_fieldnames_beh)) ~= 0        
        %number in model
        good_beh_number = find(strcmp(good_fieldnames_beh,current_beh)==1);
        %number in own dataset
        beh_num_to_change = find(strcmp(tot_fieldnames_beh,current_beh)==1);
        if good_beh_number ~= beh_num_to_change
            indices_to_change = (original_labels_newfile == beh_num_to_change) ;
            tot_annotation_output_goodtracking(indices_to_change) = good_beh_number; 
        end 
    % if current behavior has not been observed (NEED TO DECIDE HOW TO CLASSIFY THIS)                
    else
        good_fieldnames_beh = [good_fieldnames_beh; current_beh] ;
        new_beh_num = length(good_fieldnames_beh) ;
        %total is the new behavior
        beh_num_to_change = find(strcmp(tot_fieldnames_beh,current_beh)==1);
        % find from the unaltered dataset
        indices_to_change = (original_labels_newfile == beh_num_to_change) ;
        tot_annotation_output_goodtracking(indices_to_change) = new_beh_num;
    end
end