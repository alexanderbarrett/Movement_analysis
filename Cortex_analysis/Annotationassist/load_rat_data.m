function [tot_ML_features, tot_fieldnames_beh, tot_annotation_output_goodtracking, tot_annotated_frames_goodtracking] ...
    = load_rat_data(ratnames, condition_numbers_vec, annotation_numbers, mocapmasterdirectory, tot_ML_features, tot_fieldnames_beh, tot_annotation_output_goodtracking, tot_annotated_frames_goodtracking, tot_subset_of_frames_annotated)

for i = 1:length(ratnames) 
    
    ratname = char(ratnames(i)) ; 
    conditionnumbers = condition_numbers_vec(i) ; 
    annotation_number = annotation_numbers(i) ;
    
    [~,fulloutput_annotation_struct,subset_of_frames_annotated,~] = load_handannotation_files(annotation_number,mocapmasterdirectory);
    
    %% get files and the machine learning features
    descriptor_struct = get_mocap_files_table(conditionnumbers(1),ratname,mocapmasterdirectory);
    [~,mocapfilearray,~,~,~] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);      
    annotation_string = 'annotation_ii';
    ML_matobjs{1} = batch_compute_mlfeatures(mocapfilearray,ratname,0);
    
    %% load in the ML features
    [ML_features] = load_selectML(1,ML_matobjs,annotation_string);
    
    %% Convert the struct to a categorical output for the random forest, ie each behavior is a number
    subset_of_frames_annotated_full = 1:max(subset_of_frames_annotated);
    [outputvector,fieldnames_beh,~] = generate_categorical_output(fulloutput_annotation_struct,max(subset_of_frames_annotated));
    fieldnames_beh = cat(1,'null',fieldnames_beh); %there is a 0 for unannotated in the categorical vector
    
    %% remove unannotated frames
    [~,ind_unannotated] = intersect(subset_of_frames_annotated_full,find(outputvector == 0));
    subset_of_frames_annotated_full(ind_unannotated)= [];
    
    % get the outputvector- intersect the frames that have good tracking 
    [ind_for_outputvec,annotated_frames_goodtracking,~] = intersect(ML_features.framelist_true,subset_of_frames_annotated_full);
    annotation_output_goodtracking = outputvector(ind_for_outputvec);
    
    % add to aggregate data    
    % frames
    if i == 1 
        fixed_frames = annotated_frames_goodtracking ;
        fixed_subset = subset_of_frames_annotated ; 
    else
        fixed_frames = annotated_frames_goodtracking + max(tot_ML_features.framelist_true) ; 
        fixed_subset = subset_of_frames_annotated + max(tot_ML_features.framelist_true) ; 
    end
    tot_annotated_frames_goodtracking  = [tot_annotated_frames_goodtracking ; fixed_frames] ;     
    tot_subset_of_frames_annotated = [tot_subset_of_frames_annotated, fixed_subset] ;  
    
    %features 
    if numel(tot_ML_features) == 0
        tot_ML_features = ML_features ; 
    else 
        tot_ML_features = mergeStructs_AB(ML_features, tot_ML_features) ; 
    end        
    %  Labels and fieldnames_beh
    if numel(tot_fieldnames_beh) == 0
        tot_fieldnames_beh = fieldnames_beh ; 
    else  
        original_labels_newfile = annotation_output_goodtracking;
        for j = 1:length(fieldnames_beh)
            current_beh = char(fieldnames_beh(j)) ; 
            % if current behavior has already been observed
            if sum(strcmp(current_beh,tot_fieldnames_beh)) ~= 0
                good_beh_number = find(strcmp(tot_fieldnames_beh,current_beh)==1) - 1;
                beh_num_to_change = find(strcmp(fieldnames_beh,current_beh)==1) - 1;
                if good_beh_number ~= beh_num_to_change
                    indices_to_change = (original_labels_newfile == beh_num_to_change) ;
                    annotation_output_goodtracking(indices_to_change) = good_beh_number; 
                end
            % if current behavior has not been observed                
            else
                tot_fieldnames_beh = [tot_fieldnames_beh; current_beh] ;
                new_beh_num = length(tot_fieldnames_beh) - 1 ;
                beh_num_to_change = find(strcmp(fieldnames_beh,current_beh)==1) - 1;
                indices_to_change = (original_labels_newfile == beh_num_to_change) ;
                annotation_output_goodtracking(indices_to_change) = new_beh_num;                       
            end
        end
    end
    tot_annotation_output_goodtracking = [tot_annotation_output_goodtracking, annotation_output_goodtracking] ; 
end

% fix fieldnames
tot_fieldnames_beh = tot_fieldnames_beh(2:end) ; 

% remove nan
nan_indices = isnan(tot_annotation_output_goodtracking) ; 
tot_annotation_output_goodtracking = tot_annotation_output_goodtracking(~nan_indices) ; 
tot_annotated_frames_goodtracking = tot_annotated_frames_goodtracking(~nan_indices) ; 

