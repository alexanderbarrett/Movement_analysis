function make_arena_line_plots(mocapstruct)


markers_plot = [1,4,11];

fn_here = mocapstruct.markernames;

for mm = markers_plot
time_subset_line = setxor(1:size(mocapstruct.markers_preproc.(fn_here{mm}),1),mocapstruct.bad_frames_agg{mm});
figure(80+mm)
subplot(1,3,1)
line(mocapstruct.markers_preproc.(fn_here{mm})(time_subset_line,1),mocapstruct.markers_preproc.(fn_here{mm})(time_subset_line,2))
subplot(1,3,2)
line(mocapstruct.markers_preproc.(fn_here{mm})(time_subset_line,2),mocapstruct.markers_preproc.(fn_here{mm})(time_subset_line,3))
subplot(1,3,3)
line(mocapstruct.markers_preproc.(fn_here{mm})(time_subset_line,1),mocapstruct.markers_preproc.(fn_here{mm})(time_subset_line,3))



end

%get arena edges


%get axis tilts


% get lever location


%line(mocapstruct_pre.markers_preproc.SpineF(time_subset_line,1),mocapstruct_pre.markers_preproc.SpineF(time_subset_line,2))
