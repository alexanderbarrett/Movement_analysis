function [filepath_array_sorted_analog,thresh_crossing_array,num_frames_array,missing_times] = get_analogactive_files(filepath_array)

analog_active = [];
 thresh_crossing_array = zeros(1,numel(filepath_array));
  num_frames_array = zeros(1,numel(filepath_array));
  missing_times = cell(1,numel(filepath_array));

for ll = 1:numel(filepath_array)
acq = btkReadAcquisition(filepath_array{ll});

analog = btkGetAnalogs(acq);
markers = btkGetMarkers(acq);
btkCloseAcquisition(acq);

fnames = fieldnames(analog);
markernames = fieldnames(markers);

analog_threshold = 0.5;
 thresh_crossings = 0;

 
 
for mm = 1:numel(fnames)
    thresh_crossings = thresh_crossings+numel(find(analog.(fnames{mm})>analog_threshold));   
end

if (thresh_crossings >0)
    analog_active = cat(1,analog_active,ll);
    fprintf('Analog active for file %s number of crossings %f \n',filepath_array{ll},thresh_crossings);
else
        fprintf('Analog NOT active for file %s  \n',filepath_array{ll});

    
end

filepath_array_sorted_analog = filepath_array(analog_active);
thresh_crossing_array(ll) = thresh_crossings;
num_frames_array(ll) = size(markers.(markernames{1}),1);

missing_times{ll} = zeros(1,numel(markernames)); 

for lk = 1:numel(markernames)
    missing_times{ll}(lk) = numel(find(markers.(markernames{lk})(:,1)==0));
  
end


end
end