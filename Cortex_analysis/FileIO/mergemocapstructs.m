function mocap_struct_agg = mergemocapstructs(mocap_struct_agg,mocap_struct_indiv)

if ~isfield(mocap_struct_agg,'markers_preproc')
mocap_struct_agg = mocap_struct_indiv;
else
markernames = fieldnames(mocap_struct_agg.markers_preproc);
% 
agg_length = size(mocap_struct_agg.markers_preproc.(markernames{1}),1);


for mm = 1:numel(markernames)
mocap_struct_agg.markers_preproc.(markernames{mm}) = cat(1,mocap_struct_agg.markers_preproc.(markernames{mm}),...
    mocap_struct_indiv.markers_preproc.(markernames{mm}));
%mocap_struct_agg.markers_aligned_preproc.(markernames{mm}) = cat(1,mocap_struct_agg.markers_aligned_preproc.(markernames{mm}),...
 %   mocap_struct_indiv.markers_aligned_preproc.(markernames{mm}));
end

analognames = fieldnames(mocap_struct_agg.analog);
for mm = 1:numel(analognames)

mocap_struct_agg.analog.(analognames{mm}) = cat(1,mocap_struct_agg.analog.(analognames{mm}),...
    reshape(mocap_struct_indiv.analog.(analognames{mm}),[],1));
end
% 


 mocap_struct_agg.move_frames = cat(2,mocap_struct_agg.move_frames,bsxfun(@plus,mocap_struct_indiv.move_frames,agg_length));
 mocap_struct_agg.rest_frames = cat(2,mocap_struct_agg.rest_frames,bsxfun(@plus,mocap_struct_indiv.rest_frames,agg_length));
 mocap_struct_agg.move_frames_fast = cat(2,mocap_struct_agg.move_frames_fast,bsxfun(@plus,mocap_struct_indiv.move_frames_fast,agg_length));
 mocap_struct_agg.rest_frames_fast = cat(2,mocap_struct_agg.rest_frames_fast,bsxfun(@plus,mocap_struct_indiv.rest_frames_fast,agg_length));
  mocap_struct_agg.move_near_frames = cat(2,mocap_struct_agg.move_near_frames,bsxfun(@plus,mocap_struct_indiv.move_near_frames,agg_length));
 mocap_struct_agg.rest_near_frames = cat(2,mocap_struct_agg.rest_near_frames,bsxfun(@plus,mocap_struct_indiv.rest_near_frames,agg_length));
 
for mm = 1:numel(mocap_struct_agg.bad_frames_agg)
mocap_struct_agg.bad_frames_agg{mm} = cat(2,mocap_struct_agg.bad_frames_agg{mm},bsxfun(@plus,mocap_struct_indiv.bad_frames_agg{mm},agg_length));

%% remove the first and last 10 frames of files
mocap_struct_agg.bad_frames_agg{mm} = cat(2,mocap_struct_agg.bad_frames_agg{mm},agg_length+-10:10);
end

%%old messup
if (min(size(mocap_struct_indiv.resample_analog))==1)
    mocap_struct_indiv.resample_analog = reshape(mocap_struct_indiv.resample_analog,[],5);
end

mocap_struct_agg.resample_analog = cat(1,mocap_struct_agg.resample_analog,...
    reshape(mocap_struct_indiv.resample_analog,[],min(size(mocap_struct_indiv.resample_analog))));

mocap_struct_agg.lever_thresholded = cat(2,mocap_struct_agg.lever_thresholded,reshape(mocap_struct_indiv.lever_thresholded,1,[]));
mocap_struct_agg.filestartpts = cat(2,mocap_struct_agg.filestartpts, mocap_struct_indiv.filestartpts);
mocap_struct_agg.filenames = cat(2,mocap_struct_agg.filenames,mocap_struct_indiv.filenames);
mocap_struct_agg.cameradirectory = cat(1,mocap_struct_agg.cameradirectory,mocap_struct_indiv.cameradirectory);
if isfield(mocap_struct_indiv,'matched_frames_aligned')
    if numel(mocap_struct_agg.matched_frames_aligned) == 1
        mocap_struct_agg.matched_frames_aligned = mocap_struct_agg.matched_frames_aligned{1};
    end
mocap_struct_agg.matched_frames_aligned = cat(1,mocap_struct_agg.matched_frames_aligned,mocap_struct_indiv.matched_frames_aligned);
else
    mocap_struct_agg.matched_frames_aligned = cat(1,mocap_struct_agg.matched_frames_aligned,{[],[],[]});
end
%mocap_struct_agg.markernames = marker_names;
%mocap_struct_agg.fps = fps;
%mocap_struct_agg.analog_fps = analog_fps;

end
end