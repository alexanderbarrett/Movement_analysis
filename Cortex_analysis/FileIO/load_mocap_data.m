function [mocapstruct_agg] = load_mocap_data(mocapfilearray)
%load structs and aggregate
mocapstruct_agg = struct();
for mm =1:numel(mocapfilearray)
   filename_here = strrep(strrep(mocapfilearray{mm},'Generated_C3D_files\','Preprocessed\'),'.c3d','.mat');
%sort them

[mocap_struct_indiv] = load(filename_here);
if mm ==1
    mocapstruct_agg = mocap_struct_indiv.mocap_struct;
else
mocapstruct_agg = mergemocapstructs(mocapstruct_agg,mocap_struct_indiv.mocap_struct);
end
end


end