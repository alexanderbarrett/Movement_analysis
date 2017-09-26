
filepath_1 = {'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170725\Generated_C3D_files\Vicon3_recording_caffeine2_MACbatch.c3d'};
filepath_2 = {'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170725\Generated_C3D_files\Vicon3_recording_caffeine2_MACbatch_two.c3d'};
filepath_3 = {'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170725\Generated_C3D_files\Vicon3_recording_caffeine2_MACbatch_three.c3d'};
filepath_4 = {'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170725\Generated_C3D_files\Vicon3_recording_caffeine2_nolj.c3d'};


% btkCloseAcquisition(acq);
fps = 300;
analog_factor = 1;
analog_fps = fps*analog_factor;

%[markers,analog,resample_analog,lever_thresholded] = readc3d_jdm(filepath,fps,analog_fps);
[markers,analog,resample_analog,lever_thresholded] = concatenate_c3d(filepath_3,fps,analog_fps);

%% get basic information about the dataset
marker_names = fieldnames(markers);
num_markers = numel(marker_names);
% get rid of frames where markers are absent
missing_times = cell(1,numel(marker_names));
missing_times_postpreprocess = cell(1,numel(marker_names));
marker_frame_length = size(markers.(marker_names{1}),1);
markers_preproc = markers;

for ll = 1:numel(marker_names)
    missing_times{ll} = find(markers.(marker_names{ll})(:,1)==0);
    fprintf('For marker %s percentage of frames missing %e \n',marker_names{ll},...
        100.*numel(missing_times{ll})./marker_frame_length);
end

%fprintf('Percentage where two or more are gone %e \n',0)



%% other interpolation etc. 
for ll = 1:numel(marker_names)
    fprintf('starting interpolation for marker %f \n',ll)
   % for jkj = 1:3 %need to do over all x y z simul and add spike correction
    [temp,fake_frames] =markershortinterp((markers_preproc.(marker_names{ll})),30,5);
    fakedata.(marker_names{ll}){1} = fake_frames;
   % if numel(temp)
    markers_preproc.(marker_names{ll}) = temp;
  %  end
    clear fake_frames
  %  end
end

%% median filter the data to remove spikes
    fprintf('median filtering %f \n')

for ll = 1:numel(marker_names)
    markers_preproc.(marker_names{ll}) = medfilt2(markers_preproc.(marker_names{ll}),[5,1]);
end


for ll = 1:numel(marker_names)
    missing_times_postpreprocess{ll} = find(markers_preproc.(marker_names{ll})(:,1)==0);
    fprintf('For marker %s percentage of missing frames %e \n',...
        marker_names{ll},numel(missing_times_postpreprocess{ll})./marker_frame_length);
end