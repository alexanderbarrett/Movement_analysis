

[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_caff',mocapmasterdirectory);

%% either load  or preprocess from scratch
[mocapstruct] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes);


markerfields = [1,4,9,12,13,18];
fn_names = fieldnames(mocapstruct.markers_aligned_preproc);
delta = 30;
num_scales = 3;

subset = cell(1,3);
subset{1} = 1:5000000;
subset{2} = (1.5*10^6):(1.9*10^6);
subset{3} = (1.727*10^6):(1.745*10^6);

for jj = 1:num_scales
figure(400+jj)

for lk = 1:numel(markerfields)
    tracetoplot = mocapstruct.markers_aligned_preproc.(fn_names{markerfields(lk)})(:,3);
tracetoplot(mocapstruct.bad_frames_agg{markerfields(lk)}) = nan;


plot(0:1./300:(numel(subset{jj})-1)./300,delta*lk+tracetoplot(subset{jj}))
hold on
ylabel('Limb Z-Position )mm)')
xlabel('Time (s)')
end
box off
legend(fn_names(markerfields),'Location','BestOutside','Orientation','Horizontal')
end

%% animate subsets of the smallest iteration

M = animate_sbys(mocapstruct,mocapvideodirectory,subset{3})