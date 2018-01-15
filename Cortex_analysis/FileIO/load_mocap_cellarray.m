function [mocapcellout,MLfeaturescellout] = load_mocap_cellarray(ratname,numbers,mocapmasterdirectory)
mocapcellout = cell(1,numel(numbers))
MLfeaturescellout = cell(1,numel(numbers));

for ll = numbers
    cellid = find(numbers == ll);
    
descriptor_struct = get_mocap_files_table(ll,ratname);
 [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
 [mocapcellout{cellid}] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,0,0,[],mocapvideodirectory,1);
 %[mocapstruct_lesion] = preprocess_mocap_data2( mocapfilearray1,mocapfilestruct1,descriptor_struct_1,mocapfiletimes1,0,1,[],mocapvideodirectory,0);
MLfeaturescellout{cellid} = get_supervised_features(mocapcellout{cellid},mocapcellout{cellid}.modular_cluster_properties.clustering_inds_agg{2},2);
end



end