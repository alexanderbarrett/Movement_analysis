
function [marker_agg,analog_agg,resample_analog_agg,lever_thresholded_agg,filestartpts] = concatenate_c3d(filepath_array,fps,analog_fps)
marker_agg = [];
analog_agg = [];
resample_analog_agg = [];
lever_thresholded_agg = [];
filestartpts = zeros(2,numel(filepath_array));

for ll = 1:numel(filepath_array)
acq = btkReadAcquisition(filepath_array{ll});

markers = btkGetMarkers(acq);
analog = btkGetAnalogs(acq);

btkCloseAcquisition(acq);

fnames = fieldnames(analog);

%% get frame rate
fid = fopen(filepath_array{ll}, 'r');
fseek(fid,0,'cof');
frame_rate=fread(fid,1,'float32');
fclose(fid);


%% Process the analog data

w0 = 60/(analog_fps/2); %filter fundamental frequency
bw = w0/50; %Q=35
[b,a] = iirnotch(w0,bw);
if (isfield(analog,'LEVER'))
    


filtered_lever = filtfilt(b,a,analog.LEVER);

lever_threshold = 0.6;

lever_thresholded = zeros(1,numel(filtered_lever));
lever_thresholded(filtered_lever>lever_threshold) = 1;
lever_thresholded(filtered_lever<-lever_threshold) = -1;

figure(111)
plot(filtered_lever)
hold on
plot(0.6*ones(1,numel(filtered_lever)),'r')
plot(analog.LICK_SENSOR,'g')
plot(lever_thresholded,'k')
plot(analog.SPEAKER,'y')
hold off


resample_analog = resample([lever_thresholded' analog.SPEAKER analog.(fnames{3}) analog.LICK_SENSOR analog.GROUND],1,1);
resample_analog = bsxfun(@minus,resample_analog,mean(resample_analog,1));

%% supplement the length if there is a mismatch
   fn_markers = fieldnames(markers);    
marker_number = numel(numel(markers.(fn_markers{1})(:,1)));

  fn_analog = fieldnames(analog);    
analog_number = numel(analog.(fn_analog{1}));

if analog_number < marker_number
    fprintf('Mismatch for marker and analaog lengths for file number %f %s \n',ll,filepath_array{ll});
   for jj = 1:numel(fn_analog)
       analog.(fn_analog{jj}) = cat(1,analog.(fn_analog{jj}),zeros(marker_number-analog_number,1));
   end
   resample_analog = cat(1,resample_analog,zeros(marker_number-analog_number,size(resample_analog,2)));
     lever_thresholded = cat(2,lever_thresholded,zeros(size(lever_thresholded,2),marker_number-analog_number));

end

else
    
    
end


if (ll == 1)
      fn_markers = fieldnames(markers);    
      
      virtualmarkers = [];
      realmarkers = markers;
      for mm = 1:numel(fn_markers)
          if (strcmp(fn_markers{mm}(1),'V') || strcmp(fn_markers{mm},'Head_between'))
              virtualmarkers = cat(1,virtualmarkers,mm);
              realmarkers = rmfield(realmarkers,fn_markers{mm});
          end
      end
          good_markers = setxor(1:numel(fn_markers),virtualmarkers);
      
    marker_agg = realmarkers;
analog_agg = analog;
resample_analog_agg = resample_analog;
lever_thresholded_agg = lever_thresholded;

filestartpts(1,ll) = 1;
filestartpts(2,ll) = size(markers.(fn_markers{1}),1);

else
    %% concatenate analog
    fn_analog = fieldnames(analog);    
   for jj = 1:numel(fn_analog)
       analog_agg.(fn_analog{jj}) = cat(1,analog_agg.(fn_analog{jj}),analog.(fn_analog{jj}));
   end
   
      %% concatenate analog resampled
      resample_analog_agg = cat(1,resample_analog_agg,resample_analog);
%    fn_analog = fieldnames(resample_analog);    
%    for jj = 1:numel(fn_analog)
%        resample_analog_agg.(fn_analog{jj}) = cat(1,resample_analog_agg.(fn_analog{jj}),resample_analog.(fn_analog{jj}));
%    end
   
      %% concatenate markers
      filestartpts(1,ll) = size(marker_agg.(fn_markers{1}),1)+1;

   fn_markers = fieldnames(marker_agg);    
   for jj = 1:numel(fn_markers)
       marker_agg.(fn_markers{jj}) = cat(1,marker_agg.(fn_markers{jj}),markers.(fn_markers{jj}));
   end
   
   %% concatenate thresholded lever
   lever_thresholded_agg = cat(2,lever_thresholded_agg,lever_thresholded);
    
   
filestartpts(2,ll) = size(marker_agg.(fn_markers{1}),1);
   
end
end









end