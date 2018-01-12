function annotate_filelist(directory,filelist)
mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';
loadmocapfilestruct('Vicon8',mocapmasterdirectory) %this step can take an hour, potentially longer on the server
[descriptor_struct_1,mocapfilearray1,mocapfilestruct1,mocapvideodirectory,mocapfiletimes1] =  get_mocap_files('Vicon8','Vicon8_caff',mocapmasterdirectory);

% load existing random forest model for this dataset/aggregated annotation
% for this condition
annotation_xls_in = load('Y:\Jesse\Data\Annotation\annotation_categories.xls');
annotation_model_in = load('Y:\Jesse\Data\Annotation\annotation_categories.xls');

% loop over all mocap files
for mm = 1:numel(mocapfilearray1)
    %get mocap struct without saving
    [mocapstruct_in] = preprocess_mocap_data_2( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,0,0);

    % get ML features for the prediction
        ML_features = get_supervised_features(mocapstruct_caff,mocapstruct_caff.modular_cluster_properties.clustering_inds_agg{2},2);

    %call annotation with fixed output name and fixed model for condition
     MCC(annotation_xls_in,mocapstruct_in,ml_features_in,model_in); 
    %or run auto annotation

end
    


end