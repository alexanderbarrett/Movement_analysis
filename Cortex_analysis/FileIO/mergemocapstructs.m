function mocap_struct_agg = mergemocapstructs(mocap_struct_agg,mocap_struct_indiv)

markernames = fieldnames(mocap_struct_agg.markers_preproc);

for mm = 1:numel(markernames)
mocap_struct_agg.markers_preproc.(markernames{mm}) = cat(1,mocap_struct_agg.markers_preproc.(markernames{mm}),...
    mocap_struct_indiv.markers_preproc.(markernames{mm}));
mocap_struct_agg.markers_aligned_preproc.(markernames{mm}) = cat(1,mocap_struct_agg.markers_aligned_preproc.(markernames{mm}),...
    mocap_struct_indiv.markers_aligned_preproc.(markernames{mm}));
end

analognames = fieldnames(mocap_struct_agg.analog);
for mm = 1:numel(analognames)

mocap_struct_agg.analog.(analognames{mm}) = cat(1,mocap_struct_agg.analog.(analognames{mm}),...
    mocap_struct_indiv.analog.(analognames{mm}));
end

mocap_struct_agg.move_frames = cat(2,mocap_struct_agg.move_frames,mocap_struct_indiv.move_frames);
mocap_struct_agg.rest_frames = cat(2,mocap_struct_agg.rest_frames,mocap_struct_indiv.rest_frames);

for mm = 1:numel(mocap_struct_agg.bad_frames_agg)
mocap_struct_agg.bad_frames_agg{mm} = cat(2,mocap_struct_agg.bad_frames_agg{mm},mocap_struct_indiv.bad_frames_agg{mm});
end

mocap_struct_agg.fraction_missing = cat(2,mocap_struct_agg.fraction_missing,mocap_struct_indiv.fraction_missing);
mocap_struct_agg.resample_analog = cat(1,mocap_struct_agg.resample_analog,mocap_struct_indiv.resample_analog);
mocap_struct_agg.lever_thresholded = cat(2,mocap_struct_agg.lever_thresholded,mocap_struct_indiv.lever_thresholded);
mocap_struct_agg.filestartpts = cat(2,mocap_struct_agg.filestartpts, mocap_struct_indiv.filestartpts);
%mocap_struct_agg.markernames = marker_names;
%mocap_struct_agg.fps = fps;
%mocap_struct_agg.analog_fps = analog_fps;


end