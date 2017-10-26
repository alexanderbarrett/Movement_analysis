

[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_caff',mocapmasterdirectory);

%% either load  or preprocess from scratch
[mocapstruct] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes);
savedirectory = strcat(mocapmasterdirectory,'Plots_timescales_examples',filesep);
mkdir(savedirectory);

markerfields = [1,4,9,12,13,18];
fn_names = fieldnames(mocapstruct.markers_aligned_preproc);
%delta = 30;
num_scales = 3;

subset = cell(1,3);
subset{1} = 1:5000000;
subset{2} = (1.5*10^6):(1.9*10^6);
subset{3} = (1.727*10^6):(1.745*10^6);

subset{4} = (317*300):(328*300);
subset{4} = setxor(subset{4},intersect(cat(2,mocapstruct.bad_frames_agg{5},mocapstruct.bad_frames_agg{4}),subset{4}));
subset{5} = (1):(5000);


for jj = 4;%1:num_scales
figure(400+jj)

for lk = 1:numel(markerfields)
    tracetoplot = mocapstruct.markers_aligned_preproc.(fn_names{markerfields(lk)})(:,3);
tracetoplot(mocapstruct.bad_frames_agg{markerfields(lk)}) = nan;


plot(0:1./300:(numel(subset{jj})-1)./300,tracetoplot(subset{jj}))
hold on
ylabel('Limb Z-Position relative to Spine)mm)')
xlabel('Time (s)')
end
box off
legend(fn_names(markerfields)','Location','BestOutside','Orientation','Horizontal')
end
print('-depsc',strcat(savedirectory,'example_timeseries.eps'))
print('-dpng',strcat(savedirectory,'example_timeseries.png'))

%% animate subsets of the smallest iteration

v = VideoWriter(strcat(savedirectory,'sbys_movie6'),'MPEG-4');
                open(v)              
M = animate_sbys(mocapstruct,mocapvideodirectory,subset{4} )
   writeVideo(v,M)
                                    close(v)
                                    
                                    
                                    
                                    
                                    %% now get very long example
                                    
[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_caff',mocapmasterdirectory);

%% either load  or preprocess from scratch
[mocapstruct] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes);
savedirectory = strcat(mocapmasterdirectory,'Plots_timescales_examples',filesep);

                                    Vicon8_prelesion_long
                                    