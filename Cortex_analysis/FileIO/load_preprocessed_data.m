function [mocap_struct_agg] = load_preprocessed_data(filepath_array_sorted)
mocap_struct_agg = [];
for mm =1:numel(filepath_array_sorted)
    %get preproc filepath
     [save_dir,filename_here,~] = fileparts(filepath_array_sorted{mm});

    preproc_save_directory = strrep(save_dir,'Generated_C3D_files','Preprocessed\');
save_filename= strcat(preproc_save_directory,filename_here,'.mat');

    %load from disk 
  load(save_filename);
    
mocap_struct_agg = mergemocapstructs(mocap_struct_agg,mocap_struct);

end


end