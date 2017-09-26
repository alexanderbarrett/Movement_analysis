
function [marker_agg,analog_agg,resample_analog_agg,lever_thresholded_agg,fakedata,badframes,bad_frames_permarker] = concatenate_c3d(filepath_array,fps,analog_fps)
marker_agg = [];
analog_agg = [];
resample_analog_agg = [];
lever_thresholded_agg = [];

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
    analog = [];
    resample_analog = [];
    lever_thresholded = [];
end


if (numel(marker_agg)==0)
    marker_agg = markers;
else

  fn_markers = fieldnames(marker_agg);    
   for jj = 1:numel(fn_markers)
       marker_agg.(fn_markers{jj}) = cat(1,marker_agg.(fn_markers{jj}),markers.(fn_markers{jj}));
   end
   
   size_markers = size(marker_agg.(fn_markers{1}),1);
end

if (numel(analog_agg) == 0)
analog_agg = analog;
resample_analog_agg = resample_analog;
lever_thresholded_agg = lever_thresholded;

else
    %% concatenate analog
    
      %% concatenate markers
   fn_markers = fieldnames(marker_agg);    

   size_markers = size(marker_agg.(fn_markers{1}),1);
    
   
    fn_analog = fieldnames(analog_agg);    
   for jj = 1:numel(fn_analog)
       if (isfield(analog,fn_analog{jj}))
       analog_agg.(fn_analog{jj}) = cat(1,analog_agg.(fn_analog{jj}),analog.(fn_analog{jj}));
       else
            analog_agg.(fn_analog{jj}) = cat(1,analog_agg.(fn_analog{jj}),zeros( size_markers,1));
       end
   end
   
      %% concatenate analog resampled
      if numel(resample_analog)
      resample_analog_agg = cat(1,resample_analog_agg,resample_analog);
      else
            resample_analog_agg = cat(1,resample_analog_agg,zeros( size_markers,5));
    
      end
%    fn_analog = fieldnames(resample_analog);    
%    for jj = 1:numel(fn_analog)
%        resample_analog_agg.(fn_analog{jj}) = cat(1,resample_analog_agg.(fn_analog{jj}),resample_analog.(fn_analog{jj}));
%    end
   
    
   %% concatenate thresholded lever
    if numel(lever_thresholded_agg)
   lever_thresholded_agg = cat(2,lever_thresholded_agg,lever_thresholded);
    else
           lever_thresholded_agg = cat(2,lever_thresholded_agg,zeros( 1,size_markers));

    end
    
end
end


if numel(marker_agg)

fakedata = struct();
for fn = fieldnames(marker_agg)'
   fakedata.(fn{1}) = cell(1,1);
end

marker_names = fieldnames(marker_agg);
%% other interpolation etc. 
for ll = 1:numel(marker_names)
    fprintf('starting interpolation for markers, over 60 frames (200 ms) maximum %f \n',ll)
   % for jkj = 1:3 %need to do over all x y z simul and add spike correction
    [temp,fake_frames] =markershortinterp((marker_agg.(marker_names{ll})),60,5);
    fakedata.(marker_names{ll}){1} = fake_frames;
   % if numel(temp)
    marker_agg.(marker_names{ll}) = temp;
  %  end
    clear fake_frames
  %  end
end

%% median filter the data to remove spikes
    fprintf('median filtering over 3 frames %f \n')

for ll = 1:numel(marker_names)
    marker_agg.(marker_names{ll}) = medfilt2(marker_agg.(marker_names{ll}),[3,1]);
end

%remove null fields
fprintf('Removing null fields \n')

for ll = 1:numel(marker_names)
    if sum(   marker_agg.(marker_names{ll})) == 0
   marker_agg =  rmfield(marker_agg,marker_names{ll});
    end
end
marker_names = fieldnames(marker_agg);

%% now remove and identify bad frames
fprintf('Removing bad frames (will create a big matrix, watchout!!) \n')
bad_threshold = 0;
bad_frames_permarker = cell(1,numel(marker_names));
 agg_features = cell2mat(struct2cell(marker_agg')');
 
 bad_frames = find(agg_features == 0);
 [bad_inds,~] = ind2sub(size(agg_features),bad_frames);
 bad_inds = unique(bad_inds);
 badframes = bad_inds;
 
  fprintf('nubmer of inds %f number of bad frames  %f percentage %f \n', size(agg_features,1),numel(bad_inds),numel(bad_inds)./size(agg_features,1))

 
 for ll = 1:numel(marker_names)
 bad_frames_permarker{ll} = find(sum(marker_agg.(marker_names{ll}),2) == 0);
 fprintf('FOR MARKER %f nubmer of inds %f number of bad frames  %f percentage %f \n', ll,size(agg_features,1),numel(bad_frames_permarker{ll}),...
     numel(bad_frames_permarker{ll})./size(agg_features,1))
 end
 
clear agg_features
else
    fakedata = [];
    badframes = [];
     bad_frames_permarker = [];
end


end