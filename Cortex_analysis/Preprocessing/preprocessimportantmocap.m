function preprocessimportantmocap(ratname,indexes)
%1,5,6,10 for jdm 25
%1,5,6,9 for vicon8
%1,5,6,9 for jdm32

mocapmasterdirectory = 'Y:\Jesse\Data\Motionanalysis_captures\';
for mm = indexes
    fprintf('starting preprocessing for index %f \n',mm)
descriptor_struct = get_mocap_files_table(mm,ratname);
if numel(descriptor_struct)
 [~,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] = get_mocap_files_shortened(descriptor_struct,mocapmasterdirectory);
 if numel(mocapfilearray)
 [mocapstruct_v8pl2] = preprocess_mocap_data_2(mocapfilearray,mocapfilestruct,descriptor_struct,mocapfiletimes,1,1,[],mocapvideodirectory,0);
 end
end
end
end