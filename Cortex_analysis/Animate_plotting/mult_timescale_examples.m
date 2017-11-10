mocapmasterdirectory = '\\140.247.178.37\Jesse\Motionanalysis_captures\';
savedirectory = strcat(mocapmasterdirectory,'Plots',filesep);
mkdir(savedirectory);

overwrite=0;
overwrite_macro =0;
[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_amph',mocapmasterdirectory);

%% either load  or preprocess from scratch
[mocapstruct] = preprocess_mocap_data(mocapfilearray,mocapfilestruct,descriptor_struct_1,mocapfiletimes,overwrite,overwrite_macro);
%mocapstruct = mocapstruct.mocap_struct;
savedirectory = strcat(mocapmasterdirectory,'Plots_timescales_big',filesep);
mkdir(savedirectory);

markerfields = [1,4,9,12,13,18];
fn_names = fieldnames(mocapstruct.markers_aligned_preproc);
%delta = 30;
num_scales = 3;

subset = cell(1,3);
subset{1} = 1:size( mocapstruct.aligned_mean_position,1);
subset{2} = (1.5*10^6):(1.9*10^6);
subset{3} = (1.727*10^6):(1.745*10^6);

subset{4} = (317*300):(328*300);
subset{4} = setxor(subset{4},intersect(cat(2,mocapstruct.bad_frames_agg{5},mocapstruct.bad_frames_agg{4}),subset{4}));
subset{5} = (1):(5000);
subset{5} = 17*300*60:18*300*60;

for jj = 1;%1:num_scales
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
print('-depsc',strcat(savedirectory,'example_timeseries',num2str(jj),'.eps'))
print('-dpng',strcat(savedirectory,'example_timeseries',num2str(jj),'.png'))
end

for jj = 1;%1:num_scales
figure(500+jj)

    tracetoplot =  mocapstruct.aligned_mean_position(:,:);%mocapstruct.markers_aligned_preproc.(fn_names{markerfields(lk)})(:,3);
    params.fps = 300;
    markers_here = struct('fieldnamehere',[]);
    markers_here.fieldnamehere = tracetoplot;
    tracetoplot =hipass_clip(markers_here,mocapstruct.bad_frames_agg{markerfields(5)},params);
%align_hands_elbows  tracetoplot(mocapstruct.bad_frames_agg{markerfields(5)}) = 0;


[~,fr_temp,time_clustering,pc_spectrograms_temp] = spectrogram(double(tracetoplot.fieldnamehere((1:3:end),3)),mocapstruct.fps*2,...
     mocapstruct.fps,0:0.25:30,mocapstruct.fps./3);
imagesc(time_clustering./(60*60),fr_temp,log10((pc_spectrograms_temp)));
%plot(0:1./300:(numel(subset{jj})-1)./300,tracetoplot(subset{jj}))
%hold on
 caxis([-9 0])
h = colorbar;
ylabel(h, 'Log10 Spectral Power')
xlim([0 48])
ylabel('Frequency (Hz)')
xlabel('Time (hrs)')

box off
legend(fn_names(markerfields)','Location','BestOutside','Orientation','Horizontal')
print('-depsc',strcat(savedirectory,'example_positionspectrogrambig',num2str(jj),'.eps'))
print('-dpng',strcat(savedirectory,'example_positionspectrogrambig',num2str(jj),'.png'))
end







%% animate subsets of the smallest iteration

v = VideoWriter(strcat(savedirectory,'sbys_movie6'),'MPEG-4');
                open(v)              
M = animate_sbys(mocapstruct,mocapvideodirectory,subset{4} )
   writeVideo(v,M)
                                    close(v)
                                    
                                    
                                    
                                    
 %% now get very long example
                                    
[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_lesion_long',mocapmasterdirectory);
%[descriptor_struct_1,mocapfilearray,mocapfilestruct,mocapvideodirectory,mocapfiletimes] =  get_mocap_files('Vicon8','Vicon8_videotest',mocapmasterdirectory);

%% either load  or preprocess from scratch
overwrite=0;
overwrite_macro = 0;
[mocapstruct] = preprocess_mocap_data(mocapfilearray(1:50),mocapfilestruct,descriptor_struct_1,mocapfiletimes,overwrite,overwrite_macro);
savedirectory = strcat(mocapmasterdirectory,'Plots_timescales_examples_long',filesep);

                                    Vicon8_prelesion_long
                                    